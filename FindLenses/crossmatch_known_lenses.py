# Python script to compare master lensing database to DESI observations
# compares the coordinates

# imports
import sys, time
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

def main():
    start = time.time()

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--specprod', help='DESI spectroscopic production to pull from', type=str, default='iron')
    p.add_argument('--searchRadius', help='Search radius in arcseconds, default is 1', type=int, default=1)
    args = p.parse_args()
    
    # get master lens coordinates as sky coords to compare with desi observations
    masterlens = pd.read_csv('masterlens.tsv', sep='\t', skiprows=1)
    masterlensNames = masterlens.iloc[:,0].to_numpy() # get list of source names from master lens database

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

    # write out file I need for fastspecfit
    fastspecInput = desi[['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX']][desiIdx]
    fastspecInput.write('fastspec-input.fits', overwrite=True)

    # write out desi target matches
    desi[desiIdx].write(f'{args.specprod}-matches.fits', overwrite=True)

    # fix and write out masterlens matches
    masterlens.columns = masterlens.columns.str.strip('#').str.strip() # get rid of spaces
    masterlens.columns = [colname[1:-1] for colname in masterlens.columns] # get rid of extra quotes
    masterlens.columns = masterlens.columns.str.upper() # make upper case column names
    Table.from_pandas(masterlens.iloc[masterlensIdx]).write('masterlens-matches.fits', overwrite=True)

    print(f'{len(desiIdx)} matches found in {time.time()-start:.2f} seconds from the {args.specprod} production')
    
if __name__ == '__main__':
    sys.exit(main())
