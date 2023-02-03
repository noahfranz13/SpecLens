'''
SpecLens Class to hold code to separate and analyze spectroscopic
lenses in DESI spectroscopy

Author: Noah Franz
'''

import os
import numpy as np
from astropy.table import Table, vstack
import fitsio

from desispec.spectra import stack as spectraStack

class SpecLens():

    def __init__(self, filepath, outdir=os.getcwd()):
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

        self.infile = Table(fitsio.read(self.filepath))[:10] # FIX ME: only grabbing first 10 lines for testing

        # initiate some instance variables for later
        self.zbestfile = None
        self.coadds = []
        
    def _preprocess(self, specprod='fuji'):
        """Can't use the Docker container for pre-processing because we do not have
        redrock:

        fastspecfit is aliased to `source $IMPY_DIR/bin/fastspecfit-env-nersc'

        """
        from redrock.external.desi import write_zbest
        from desispec.io import write_spectra, read_spectra
        
        # extract info from z catalog
        zbests = []
        fibermaps = []
        expfibermaps = []
        tsnr2s = []
        for row in self.infile:

            targetid = row['TARGETID']
            survey = row['SURVEY'] 
            program = row['PROGRAM']
            hpx = row['HEALPIX']
            
            datadir = f'/global/cfs/cdirs/desi/spectro/redux/fuji/healpix/{survey}/{program}/{str(hpx)[:-2]}/{hpx}'

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
            
            spechdr = fitsio.read_header(coaddfile)
            spec = read_spectra(coaddfile).select(targets=targetid)
            assert(np.all(spec.fibermap['TARGETID'] == targetid))

            # update the headers so things work with fastspecfit
            spechdr['SPGRP'] = 'healpix'
            spechdr['SPGRPVAL'] = 0
            spechdr['SURVEY'] = 'speclens' #survey 
            spechdr['PROGRAM'] = 'speclens' #program
            spechdr['SPECPROD'] = specprod
            spec.meta = spechdr
            
            self.coadds.append(spec)
            
            zbests.append(zbest)
            fibermaps.append(fibermap)
            expfibermaps.append(expfibermap)
            tsnr2s.append(tsnr2)

        # update the headers so things work with fastspecfit
        redhdr['SPGRP'] = 'healpix'
        redhdr['SPGRPVAL'] = 0
        redhdr['SURVEY'] = 'speclens' #survey
        redhdr['PROGRAM'] = 'speclens' #program
        redhdr['SPECPROD'] = specprod

        # write out zbests
        self.zbestfile = os.path.join(self.outdir, 'redrock-preprocess.fits')
        archetype_version = None
        template_version = {redhdr['TEMNAM{:02d}'.format(nn)]: redhdr['TEMVER{:02d}'.format(nn)] for nn in np.arange(10)}
        print(f'Writing {self.zbestfile}')
        zbest = vstack(zbests, join_type='exact')
        
        # for now just take the first of each of these
        fibermap = fibermaps[0]#vstack(fibermaps, join_type='exact')
        expfibermap = expfibermaps[0] #vstack(expfibermaps, join_type='exact')
        tsnr2 = tsnr2s[0] #vstack(tsnr2s, join_type='exact')

        write_zbest(self.zbestfile, zbest, fibermap, expfibermap, tsnr2,
                    template_version, archetype_version, spec_header=redhdr)
        
    def modelLens(self, makeqa=False, mp=1):
        '''
        Model the lens (hopefully) using fastspecfit

        Note that if the source is not strong enough 
        this may pick up the entire spectrum, not
        just the lens
        '''

        # imports
        from fastspecfit.mpi import plan
        from fastspecfit.continuum import ContinuumTools
        from fastspecfit.emlines import EMLineFit
        from fastspecfit.io import DESISpectra, write_fastspecfit, read_fastspecfit
        from fastspecfit.fastspecfit import _fastspec_one, fastspec_one, _desiqa_one, desiqa_one
        
        fastfitfile = os.path.join(self.outdir, 'fastspec-speclens.fits')
        
        targetids = self.infile['TARGETID'] 
        
        self._preprocess() # preprocess the data

        sample = Table(fitsio.read(self.zbestfile))
        sample = sample[np.isin(sample['TARGETID'], targetids)]

        Spec = DESISpectra()
        CFit = ContinuumTools()
        EMFit = EMLineFit()

        if makeqa:
            print('Making QA Plots...')

            qadir = os.path.join(self.outdir, 'qa')
            if not os.path.isdir(qadir):
                os.makedirs(qadir, exist_ok=True)
            
            fastfit, metadata, coadd_type, _ = read_fastspecfit(fastfitfile)
            targetids = metadata['TARGETID'].data

            Spec.select(redrockfiles=self.zbestfile, targetids=targetids, use_quasarnet=False)
            data = Spec.read_and_unpack(CFit, fastphot=False, synthphot=True, remember_coadd=True)
        
            indx = np.arange(len(data))
            qaargs = [(CFit, EMFit, data[igal], fastfit[indx[igal]], metadata[indx[igal]],
                       coadd_type, False, qadir, None) for igal in np.arange(len(indx))]                

            if mp > 1:
                import multiprocessing
                with multiprocessing.Pool(mp) as P:
                    P.map(_desiqa_one, qaargs)
            else:
                [desiqa_one(*_qaargs) for _qaargs in qaargs]
        else:
            Spec.select(redrockfiles=self.zbestfile, targetids=targetids, zmin=-0.1,
                        use_quasarnet=False, ntargets=len(self.infile))
            data = Spec.read_and_unpack(CFit, fastphot=False, synthphot=True, remember_coadd=True)

            out, meta = Spec.init_output(CFit=CFit, EMFit=EMFit, fastphot=False)

            # Fit in parallel
            t0 = time.time()
            fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit, EMFit, False, False) # verbose and broadlinefit
                       for iobj in np.arange(Spec.ntargets)]
            if mp > 1:
                import multiprocessing
                with multiprocessing.Pool(mp) as P:
                    _out = P.map(_fastspec_one, fitargs)
            else:
                _out = [fastspec_one(*_fitargs) for _fitargs in fitargs]
            _out = list(zip(*_out))
            out = Table(np.hstack(_out[0]))
            meta = Table(np.hstack(_out[1]))

            try:
                # need to vstack to preserve the wavelength metadata 
                modelspec = vstack(_out[2], metadata_conflicts='error')
            except:
                errmsg = 'Metadata conflict when stacking model spectra.'
                log.critical(errmsg)
                raise ValueError(errmsg)

            log.info('Fitting everything took: {:.2f} sec'.format(time.time()-t0))

            # Write out.
            write_fastspecfit(out, meta, outfile=fastfitfile, modelspectra=modelspec, specprod=Spec.specprod,
                          coadd_type=Spec.coadd_type, fastphot=False)


        def modelSource(self):
            '''
            Model the source galaxy (the one being lensed) by
            subtracting the lens model from the original spectra

            Note that if the features of the source spectra are not 
            strong enough then there will be a featureless 
            spectrum leftoover.
            '''
            pass
