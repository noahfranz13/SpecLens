'''
SpecLens Class to hold code to separate and analyze spectroscopic
lenses in DESI spectroscopy

Author: Noah Franz
'''

import os, pdb
from subprocess import run

import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import fitsio

from desispec.io import write_spectra, read_spectra        

class SpecLens():

    def __init__(self, filepath, outdir=None, specprod='iron', overwrite=True, mp=1):
        '''
        filepath [string] : path to fits file with the DESI Spectroscopic Lens
                            observations. Must have the at least the following 
                            columns:
                            1) TARGETID
                            2) SURVEY
                            3) PROGRAM
                            4) HEALPIX
        outdir [string]   : output directory for files written by this class, 
                            default is None and sets to the directory of the
                            input filepath.
        specprod [string] : spectroscopic production to run on
        overwrite [bool]  : if True, overwrites existing files, default is True
        mp [int]          : number of cores to use with multiprocessing, default is 1
        '''
        
        self.filepath = filepath
        if outdir is None:
            self.outdir = os.path.split(self.filepath)[0]
        else:
            self.outdir = outdir
        self.specprod = specprod

        self.ntest = 3

        self.infile = Table(fitsio.read(self.filepath))[:self.ntest] # FIX ME: only grabbing first 10 lines for testing

        # initiate some paths for later
        self.zbestfile = os.path.join(self.outdir, 'redrock-lenses.fits')
        self.coaddfile = os.path.join(self.outdir, 'coadd-lenses.fits')
        self.fastspecfile = os.path.join(self.outdir, 'fastspec-lenses.fits')
        self.sourcespecfile = os.path.join(self.outdir, 'coadd-source.fits')
        self.sourcezbestfile = os.path.join(self.outdir, 'redrock-source.fits')
        self.QAdir = os.path.join(self.outdir, 'QA')
        
        if overwrite or (not os.path.isfile(self.zbestfile) and not os.path.isfile(self.coaddfile)):
            self.spectra = None
            self._preprocess() # preprocess the data
        else:
            self.spectra = read_spectra(self.coaddfile)
            
        # save stacked coadd spectra
        self.sourceSpec = None

        # set up multiprocessing
        self.mp = mp
        
    def _preprocess(self):
        """Can't use the Docker container for pre-processing because we do not have
        redrock:

        fastspecfit is aliased to `source $IMPY_DIR/bin/fastspecfit-env-nersc'

        """
        from redrock.external.desi import write_zbest
        from desispec.spectra import stack

        print('\n\nStarting preprocessing...')
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
        self.spectra = stack(coadds)
        print(self.spectra)
        write_spectra(self.coaddfile, self.spectra)
        
    def modelLens(self):
        '''
        Model the lens (hopefully) using fastspecfit

        Note that if the source is not strong enough 
        this may pick up the entire spectrum, not
        just the lens
        '''
        
        print('\n\nModelling the lens with fastspecfit...')
        
        fastfitfile = self.fastspecfile

        # add arguments to a list for subprocess
        arg = ['fastspec', '--mp', f'{self.mp}',
                           '--outfile', f'{fastfitfile}',
                           '--specproddir', f'{self.specprod}',
                           f'{self.zbestfile}']
        run(arg)
        
    def modelSource(self, overwrite=False):
        '''
        Model the source galaxy (the one being lensed) by
        subtracting the lens model from the original spectra

        Note that if the features of the source spectra are not 
        strong enough then there will be a featureless 
        spectrum leftover.
        '''
        print('\n\nModeling the source...') 
        self.subtract()

        from redrock.external.desi import rrdesi
            
        if not overwrite and os.path.exists(self.sourcezbestfile):
            print('Overwrite set to False and redrock file exists so skipping...')
        else:
            rrdesi(options=['--mp', str(self.mp), '-o', self.sourcezbestfile, '-i', self.sourcespecfile])
        
    def subtract(self):
        '''
        Subtracts the lens model from the combined spectra
        '''
        from copy import deepcopy
        from desispec.io import write_spectra
        
        lensWave, lensFlux, lensMeta = self.getFastspecModel()
        self.prepSpectra(lensMeta)

        combFlux = self.spectra.flux['brz']
        assert np.all(np.isclose(lensWave, self.spectra.wave['brz']))

        # subtract the fluxes
        print(len(combFlux), len(lensFlux))
        print(combFlux, lensFlux)
        sourceFlux = combFlux - lensFlux
        self.sourceSpec = deepcopy(self.spectra)

        self.sourceSpec.flux['brz'] = sourceFlux
        #pdb.set_trace()
        
        # write out spectra object
        write_spectra(self.sourcespecfile, self.sourceSpec)

    def getFastspecModel(self):
        '''
        Reads in fastspec model and preps it for subtraction
        '''
        
        # read in the metadata
        modelmeta = Table(fitsio.read(self.fastspecfile, 'METADATA'))
        
        # read in fastspecfit model
        models, hdr = fitsio.read(self.fastspecfile, 'MODELS', header=True)
        modelwave = hdr['CRVAL1'] + np.arange(hdr['NAXIS1']) * hdr['CDELT1']
    
        modelflux = [np.sum(m, axis=0).flatten() for m in models] # THIS LINE IS BUGGED WITH MULTIPLE OBJECTS!
        
        return modelwave, modelflux, modelmeta

    def prepSpectra(self, fastmeta):
        '''
        Preprocesses the Spectra and collapse b, r, z bands
        into one for modelling

        fastmeta : fastspecfit model meta data
        '''
        from desispec.coaddition import coadd_cameras
        from desiutil.dust import dust_transmission
        
        coaddSpec = coadd_cameras(self.spectra)
        bands = coaddSpec.bands[0]
        # milky way dust transmission
        mwSpec = np.array([dust_transmission(coaddSpec.wave[bands], row['EBV']) for row in fastmeta])

        # make the correction and save update spectra in instance variable
        coaddSpec.flux = {'brz' : coaddSpec.flux[bands] / mwSpec}
        self.spectra = coaddSpec

    def separateLens(self):
        '''
        Uses the rest of the methods in this class to separate the spectra
        (Basically just does everything that we could need)
        '''

        self.modelLens()
        self.modelSource()

    def generateQA(self):
        '''
        Wrapper function to generate quality analysis (QA) plots for the 
        fastspecfit run during the lens modelling step of the code. If 
        fastspecfit has not been run it runs it first and then generates 
        the plots.
        '''
        if not os.path.isfile(self.fastspecfile):
            # run fastspecfit if the file hasn't been generated
            self.modelLens()

        if not os.path.exists(self.QAdir):
            os.mkdir(self.QAdir)
        
        cmd = ['fastspecfit-qa', self.fastspecfile,
               '-o', self.QAdir,
               '--redrockfiles', self.zbestfile,
               '--mp', str(self.mp)]
        run(cmd)
