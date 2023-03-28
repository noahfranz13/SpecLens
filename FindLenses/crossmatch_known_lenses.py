# Python script to compare master lensing database to DESI observations
# compares the coordinates
# imports
import os, sys, time, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord, concatenate
import astropy.units as u
from astropy.table import Table, vstack
from astropy.io import ascii
import fitsio

# useful functions
def sepHist(masterLensCoords, desiCoords):
    '''
    masterLensCoords [SkyCoord Obj] : All of the master lens coordinates
    desiCoords [SkyCoord Obj] : All of the coordinats in the DESI SpecProd Best catalog
    '''
    # calculate all of the nearest neighbors
    _, sep2d, _ = masterLensCoords.match_to_catalog_sky(desiCoords)

    # plot them up to look for a good maximum search radius
    fig, ax = plt.subplots()
    _ = ax.hist(sep2d.arcsec, bins=100, color='cornflowerblue', range=(0,5))
    ax.set_ylabel('N')
    ax.set_xlabel('Target Separation (arcsec)')

    fig.savefig('separation-hist.png', bbox_inches='tight', transparent=False)

def mask(a):
    '''
    Mask the input to drop duplicates and get it ready for fastspecfit
    '''
    
    I = np.where((a['TARGETID'] >= 0) * (a['Z'] > 1e-3) * (a['OBJTYPE'] == 'TGT') * (a['ZWARN'] & 2**9 == 0))[0] # get rid of negative targs
    a = a[I]

    print(a[['TARGETID', 'Z', 'OBJTYPE', 'ZWARN']])

    return a, I
        
def main():
    start = time.time()

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--specprod', help='DESI spectroscopic production to pull from', type=str, default='iron')
    p.add_argument('--searchRadius', help='Search radius in arcseconds, default is 1', type=int, default=1)
    p.add_argument('--outdir', default=os.getcwd(), type=str, help='output file directory')
    p.add_argument('--overwrite', action='store_true', help='Overwrite any existing output files.')
    args = p.parse_args()

    outpath = os.path.join(args.outdir, 'fastspec-input.fits')
    if not os.path.isfile(outpath) or args.overwrite:
        # get master lens coordinates as sky coords to compare with desi observations
        masterlens = pd.read_csv('masterlens.tsv', sep='\t', skiprows=1, index_col=False)
        #print(masterlens)
        #masterlensNames = masterlens[:,0].as_array() # get list of source names from master lens database
        masterlens = Table.from_pandas(masterlens)
        
        raDeg = np.array(masterlens[' "ra_coord"']).astype(float)
        decDeg = np.array(masterlens[' "dec_coord"']).astype(float)

        raDeg = raDeg[~np.isnan(raDeg)]
        decDeg = decDeg[~np.isnan(decDeg)]

        lensCoords = SkyCoord(raDeg*u.deg, decDeg*u.deg)
    
        # read in specprod catalog
        desiPath = f'/global/cfs/cdirs/desi/spectro/redux/{args.specprod}/zcatalog/zall-pix-{args.specprod}.fits'
        desi = Table(fitsio.FITS(desiPath)[1].read())
        desi = desi[desi['ZCAT_PRIMARY']] # only select the best among duplicates
        
        # convert the specprod ra and dec to SkyCoords
        desiCoords = SkyCoord(desi['TARGET_RA']*u.deg, desi['TARGET_DEC']*u.deg)

        # plot up histogram of separations
        sepHist(lensCoords, desiCoords)
    
        # perform the crossmatch
        searchSep = args.searchRadius*u.arcsec # based on the histogram, 1 arcsecond seems reasonable
        desiIdx, masterlensIdx, sep, _ = lensCoords.search_around_sky(desiCoords, searchSep)
            
        # drop duplicates
        desiout = desi[desiIdx]
        masterlensout = masterlens[masterlensIdx]
        _, I = np.unique(desiout['TARGETID'], return_index=True)
        desiout = desiout[I]
        masterlensout = masterlensout[I]

        desiout, I2 = mask(desiout)
        masterlensout = masterlensout[I2]

        # write out the separation array
        sep = np.array(sep[I][I2])
        np.savetxt(f'{args.specprod}-matches-separation.txt', sep)
        
        # write out file I need for fastspecfit
        fastspecInput = desiout[['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX']]
        fastspecInput.write(outpath, overwrite=True)

        # write out desi target matches
        desiout.write(os.path.join(args.outdir, f'{args.specprod}-matches.fits'), overwrite=True)

        # fix and write out masterlens matches
        tout = Table(names=[n.strip('#').strip()[1:-1].upper() for n in masterlensout.columns], dtype=masterlensout.dtype)
        for row in masterlensout:
            tout.add_row(row)
        tout.write(os.path.join(args.outdir, 'masterlens-matches.fits'), overwrite=True)

        print(f'{len(desiout)} matches found in {time.time()-start:.2f} seconds from the {args.specprod} production')
    else:
        print('Output file found so not running crossmatch!')
        
if __name__ == '__main__':
    sys.exit(main())
