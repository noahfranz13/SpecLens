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
def sepHist(masterLensCoords, fujiCoords):
    '''
    masterLensCoords [SkyCoord Obj] : All of the master lens coordinates
    fujiCoords [SkyCoord Obj] : All of the coordinats in the Fuji Best catalog
    '''
    # calculate all of the nearest neighbors
    _, sep2d, _ = masterLensCoords.match_to_catalog_sky(fujiCoords)

    # plot them up to look for a good maximum search radius
    fig, ax = plt.subplots()
    _ = ax.hist(sep2d.arcsec, bins=100, color='cornflowerblue', range=(0,5))
    ax.set_ylabel('N')
    ax.set_xlabel('Target Separation (arcsec)')

    fig.savefig('separation-hist.png', bbox_inches='tight', transparent=False)

def main():
    start = time.time()
    
    # get master lens coordinates as sky coords to compare with desi observations
    masterlens = pd.read_csv('masterlens.tsv', sep='\t', skiprows=1)
    masterlensNames = masterlens.iloc[:,0].to_numpy() # get list of source names from master lens database

    raDeg = np.array(masterlens[' "ra_coord"']).astype(float)
    decDeg = np.array(masterlens[' "dec_coord"']).astype(float)

    raDeg = raDeg[~np.isnan(raDeg)]
    decDeg = decDeg[~np.isnan(decDeg)]

    lensCoords = SkyCoord(raDeg*u.deg, decDeg*u.deg)
    
    # read in fuji catalog
    fujiPath = '/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zall-pix-fuji.fits'
    fuji = Table(fitsio.FITS(fujiPath)[1].read())
    fuji = fuji[fuji['ZCAT_PRIMARY']] # only select the best among duplicates
    
    # convert the fuji ra and dec to SkyCoords
    fujiCoords = SkyCoord(fuji['TARGET_RA']*u.deg, fuji['TARGET_DEC']*u.deg)

    # plot up histogram of separations
    sepHist(lensCoords, fujiCoords)
    
    # perform the crossmatch
    searchSep = 1*u.arcsec # based on the histogram, 1 arcsecond seems reasonable
    fujiIdx, masterlensIdx, sep, _ = lensCoords.search_around_sky(fujiCoords, searchSep)

    # write out file I need for fastspecfit
    fastspecInput = fuji[['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX']][fujiIdx]
    fastspecInput.write('fastspec-input.fits', overwrite=True)

    # write out fuji matches
    fuji[fujiIdx].write('fuji-matches.fits', overwrite=True)

    # fix and write out masterlens matches
    masterlens.columns = masterlens.columns.str.strip('#').str.strip() # get rid of spaces
    masterlens.columns = [colname[1:-1] for colname in masterlens.columns] # get rid of extra quotes
    masterlens.columns = masterlens.columns.str.upper() # make upper case column names
    Table.from_pandas(masterlens.iloc[masterlensIdx]).write('masterlens-matches.fits', overwrite=True)

    print(f'{len(fujiIdx)} matches found in {time.time()-start:.2f} seconds')
    
if __name__ == '__main__':
    sys.exit(main())
