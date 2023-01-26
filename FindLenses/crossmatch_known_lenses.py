# Python script to compare master lensing database to DESI observations
# compares the coordinates

# imports
import sys
import pandas as pd
import numpy as np

from astropy.coordinates import SkyCoord, concatenate
import astropy.units as u
from astropy.table import Table, vstack
from astropy.io import ascii
import fitsio

# useful functions
def checkSep(desiCoords, lensCoord, desiTable, epsilon=1e-4):

    sep = np.array(lensCoord.separation(desiCoords))
    whereClose = np.where(sep < epsilon)[0]
    
    goodTarg = None
    
    if len(whereClose) != 0:
        goodTarg = desiTable[whereClose]
        
    return goodTarg

def main():
    
    # get master lens coordinates as sky coords to compare with desi observations
    masterlens = pd.read_csv('masterlens.tsv', sep='\t', skiprows=1)
    masterlensNames = masterlens.iloc[:,0].to_numpy() # get list of source names from master lens database

    raDeg = np.array(masterlens[' "ra_coord"']).astype(float)
    decDeg = np.array(masterlens[' "dec_coord"']).astype(float)

    raDeg = raDeg[~np.isnan(raDeg)]
    decDeg = decDeg[~np.isnan(decDeg)]

    lensCoords = SkyCoord(raDeg*u.deg, decDeg*u.deg)
    
    # get master list of all fuji and gaudalupe targets
    # read in fuji catalog
    fujiPath = '/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zall-pix-fuji.fits'
    fuji = Table(fitsio.read(fujiPath))

    # read in guadalupe catalog
    guadalupePath = '/global/cfs/cdirs/desi/spectro/redux/guadalupe/zcatalog/zall-pix-guadalupe.fits'
    guadalupe = Table(fitsio.read(guadalupePath))
    
    # convert the fuji and guadalupe ra and dec to SkyCoords
    fujiCoords = SkyCoord(fuji['TARGET_RA']*u.deg, fuji['TARGET_DEC']*u.deg)
    guadalupeCoords = SkyCoord(guadalupe['TARGET_RA']*u.deg, guadalupe['TARGET_DEC']*u.deg)
    
    # perform the crossmatch
    # THIS TAKES A WHILE
    
    count = 0
    coords = [fujiCoords, guadalupeCoords]
    tables = [fuji, guadalupe]
    names = ['fuji', 'guadalupe']
    goodRows = []
    for coordSet, table, filename in zip(coords, tables, names):
        
        print(f'----- {filename} -----')
        
        matches = []
        for ii, lensCoord in enumerate(lensCoords):

            match = checkSep(coordSet, lensCoord, table)
            if match is not None:
                match['OBJNAME'] = [masterlensNames[ii]]*len(match)
                matches.append(match)

                count += len(match)
                print(f'found {count}')
                
            if ii%1000 == 0:
                print(f'Index: {ii}')
        
        if len(matches) > 0:
            stack = vstack(matches)
            
            out = stack[['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'Z']]
            out.write(f'{filename}-masterlens-matches-info.fits', overwrite=True)

            ascii.write(stack, filename+'-matches.ecsv', format='ecsv', overwrite=True)
            
        print()
               
if __name__ == '__main__':
    sys.exit(main())