'''
SpecLens Class to hold code to separate and analyze spectroscopic
lenses in DESI spectroscopy

Author: Noah Franz
'''

import os, pdb
import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import fitsio

from desispec.spectra import stack as spectraStack

class SpecLens():

    def __init__(self, filepath, outdir=os.getcwd(), specprod='iron'):
        '''
        filepath [string] : path to fits file with the DESI Spectroscopic Lens
                            observations. Must have the at least the following 
                            columns:
                            1) TARGETID
                            2) SURVEY
                            3) PROGRAM
                            4) HEALPIX
        '''
        
        self.filepath = filepath
        self.outdir = outdir
        self.specprod = specprod

        self.ntest = 3
        
        self.infile = Table(fitsio.read(self.filepath))[:self.ntest] # FIX ME: only grabbing first 10 lines for testing

        # initiate some instance variables for later
        self.zbestfile = os.path.join(self.outdir, 'redrock-lenses.fits')
        self.coaddfile = os.path.join(self.outdir, 'coadd-lenses.fits')
        self.fastspecfile = os.path.join(self.outdir, 'fastspec-lenses.fits')
        
    def _preprocess(self):
        """Can't use the Docker container for pre-processing because we do not have
        redrock:

        fastspecfit is aliased to `source $IMPY_DIR/bin/fastspecfit-env-nersc'

        """
        from redrock.external.desi import write_zbest
        from desispec.io import write_spectra, read_spectra
        from desispec.spectra import stack
        
        # extract info from z catalog
        zbests = []
        fibermaps = []
        expfibermaps = []
        tsnr2s = []
        coadds = []
        for row in self.infile:

            targetid = row['TARGETID']
            survey = row['SURVEY'] 
            program = row['PROGRAM']
            hpx = row['HEALPIX']
            
            datadir = f'/global/cfs/cdirs/desi/spectro/redux/{self.specprod}/healpix/{survey}/{program}/{str(hpx)[:-2]}/{hpx}'

            # read the redrock and coadd catalog
            coaddfile = os.path.join(datadir, f'coadd-{survey}-{program}-{hpx}.fits')
            redrockfile = os.path.join(datadir, f'redrock-{survey}-{program}-{hpx}.fits')

            redhdr = fitsio.read_header(redrockfile)
            zbest = Table.read(redrockfile, 'REDSHIFTS')
            fibermap = Table.read(redrockfile, 'FIBERMAP')
            expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
            tsnr2 = Table.read(redrockfile, 'TSNR2')

            I = np.where(zbest['TARGETID'] == targetid)

            zbest = zbest[I]
            fibermap = fibermap[I]
            expfibermap = expfibermap[I]
            tsnr2 = tsnr2[I]
            assert(np.all(zbest['TARGETID'] == targetid))
            
            spechdr = fits.getheader(coaddfile)
            #spechdr = fitsio.read_header(coaddfile)
            spec = read_spectra(coaddfile).select(targets=targetid)
            assert(np.all(spec.fibermap['TARGETID'] == targetid))

            # update the headers so things work with fastspecfit
            spechdr['SPGRP'] = 'custom'
            spechdr['SPGRPVAL'] = 0
            spechdr['SURVEY'] = 'speclens' #survey 
            spechdr['PROGRAM'] = 'speclens' #program
            spechdr['SPECPROD'] = self.specprod
            spec.meta = spechdr

            coadds.append(spec)
            
            zbests.append(zbest)
            fibermaps.append(fibermap)
            expfibermaps.append(expfibermap)
            tsnr2s.append(tsnr2)

        # update the headers so things work with fastspecfit
        redhdr['SPGRP'] = 'healpix'
        redhdr['SPGRPVAL'] = 0
        redhdr['SURVEY'] = 'speclens' #survey
        redhdr['PROGRAM'] = 'speclens' #program
        redhdr['SPECPROD'] = self.specprod

        # write out zbests
        archetype_version = None
        template_version = {redhdr['TEMNAM{:02d}'.format(nn)]: redhdr['TEMVER{:02d}'.format(nn)] for nn in np.arange(self.ntest)}
        zbest = vstack(zbests, join_type='exact')
        
        # for now just take the first of each of these
        fibermap = vstack(fibermaps)#, join_type='exact')
        expfibermap = vstack(expfibermaps)#, join_type='exact')
        tsnr2 = vstack(tsnr2s)#, join_type='exact')

        print('Writing {}'.format(self.zbestfile))
        write_zbest(self.zbestfile, zbest, fibermap, expfibermap, tsnr2,
                    template_version, archetype_version, spec_header=redhdr)

        print('Writing {}'.format(self.coaddfile))
        write_spectra(self.coaddfile, stack(coadds))
        
    def modelLens(self, makeqa=False, mp=1):
        '''
        Model the lens (hopefully) using fastspecfit

        Note that if the source is not strong enough 
        this may pick up the entire spectrum, not
        just the lens
        '''

        # imports
        from subprocess import run
        
        fastfitfile = self.fastspecfile

        if not os.path.isfile(self.zbestfile) and not os.path.isfile(self.coaddfile):
            self._preprocess() # preprocess the data

        # add arguments to a list for subprocess
        arg = ['fastspec', '--mp', f'{mp}',
                           '--outfile', f'{fastfitfile}',
                           '--specproddir', f'{self.specprod}',
                           f'{self.zbestfile}']
        run(arg)
        
    def modelSource(self):
        '''
        Model the source galaxy (the one being lensed) by
        subtracting the lens model from the original spectra

        Note that if the features of the source spectra are not 
        strong enough then there will be a featureless 
        spectrum leftoover.
        '''
        pass
