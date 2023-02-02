#!/usr/bin/env python
"""
Wrapper to fit spectroscopic lenses with fastspecfit 

"""
#mport pdb # for debugging

import os, time, subprocess
import numpy as np
import fitsio
from glob import glob
from astropy.table import Table, vstack

from desiutil.log import get_logger
log = get_logger()

def preprocess(outdir, zcat, specprod='fuji'):
    """Can't use the Docker container for pre-processing because we do not have
    redrock:
   
    fastspecfit is aliased to `source $IMPY_DIR/bin/fastspecfit-env-nersc'
    
    """
    from redrock.external.desi import write_zbest
    from desispec.io import write_spectra, read_spectra
    
    # extract info from z catalog
    targetid = zcat['TARGETID']
    survey = zcat['SURVEY'][0] # just grab first value for now
    program = zcat['PROGRAM'][0] # same here

    datadir = '/global/cfs/cdirs/desi/spectro/redux/fuji/healpix/sv3/bright/259/25964'

    # read the redrock and coadd catalog
    coaddfile = os.path.join(datadir, 'coadd-sv3-bright-25964.fits')
    redrockfile = os.path.join(datadir, 'redrock-sv3-bright-25964.fits')

    outcoaddfile = os.path.join(outdir, 'coadd-sv3-bright-25964-out.fits')
    outredrockfile = os.path.join(outdir, 'redrock-sv3-bright-25964-out.fits')

    redhdr = fitsio.read_header(redrockfile)
    zbest = Table.read(redrockfile, 'REDSHIFTS')
    fibermap = Table.read(redrockfile, 'FIBERMAP')
    expfibermap = Table.read(redrockfile, 'EXP_FIBERMAP')
    tsnr2 = Table.read(redrockfile, 'TSNR2')

    I = np.hstack([np.where(tid == zbest['TARGETID'])[0] for tid in targetid])

    zbest = zbest[I]
    fibermap = fibermap[I]
    expfibermap = expfibermap[I]
    tsnr2 = tsnr2[I]
    assert(np.all(zbest['TARGETID'] == targetid))

    spechdr = fitsio.read_header(coaddfile)
    spec = read_spectra(coaddfile).select(targets=targetid)
    assert(np.all(spec.fibermap['TARGETID'] == targetid))
    
    # update the headers so things work with fastspecfit
    redhdr['SPGRP'] = 'healpix'
    redhdr['SPGRPVAL'] = 0
    redhdr['SURVEY'] = survey
    redhdr['PROGRAM'] = program
    redhdr['SPECPROD'] = specprod
    
    spechdr['SPGRP'] = 'healpix'
    spechdr['SPGRPVAL'] = 0
    spechdr['SURVEY'] = survey
    spechdr['PROGRAM'] = program
    spechdr['SPECPROD'] = specprod
    spec.meta = spechdr

    print('Writing {}'.format(outcoaddfile))
    write_spectra(outcoaddfile, spec)
    
    archetype_version = None
    template_version = {redhdr['TEMNAM{:02d}'.format(nn)]: redhdr['TEMVER{:02d}'.format(nn)] for nn in np.arange(10)}

    print('Writing {}'.format(outredrockfile))
    write_zbest(outredrockfile, zbest, fibermap, expfibermap, tsnr2,
                template_version, archetype_version, spec_header=redhdr)

    #db.set_trace()
    
def main():
    """Main wrapper on fastspec.

    """
    import argparse    
    from fastspecfit.mpi import plan
    from fastspecfit.continuum import ContinuumFit
    from fastspecfit.emlines import EMLineFit
    from fastspecfit.io import DESISpectra, write_fastspecfit, read_fastspecfit
    from fastspecfit.fastspecfit import _fastspec_one, fastspec_one, _desiqa_one, desiqa_one
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--mp', type=int, default=1, help='Number of multiprocessing processes per MPI rank or node.')
    parser.add_argument('--zcatfile', type=str, default=None, help='Path to z-catalog file with observations to process. Should have columns for survey, program, healpix, and targetid')
    
    parser.add_argument('--preprocess', action='store_true', help='Preprocess the files.')
    parser.add_argument('--makeqa', action='store_true', help='Build QA in parallel.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite any existing output files.')

    args = parser.parse_args()

    # project parameters
    outdir = '/global/cfs/cdirs/desi/users/nrfran/speclens/'
    redrockfile = os.path.join(outdir, 'redrock-sv3-bright-25964-out.fits')
    fastfitfile = os.path.join(outdir, 'fastspec-sv3-bright-25964-out.fits')
    
    qadir = os.path.join(outdir, 'qa')
    if not os.path.isdir(qadir):
        os.makedirs(qadir, exist_ok=True)
        
    # select a subset of targets
    if args.zcatfile:
        zcat = Table.read(args.zcatfile)
        targetids = zcat['TARGETID'] 
    else:
        raise IOError('Please provide a z catalog file with necessary info')

    if args.preprocess or not os.path.isfile(redrockfile):
        preprocess(outdir, zcat)

    sample = Table(fitsio.read(redrockfile))
    sample = sample[np.isin(sample['TARGETID'], targetids)]


    Spec = DESISpectra()
    CFit = ContinuumFit()
    EMFit = EMLineFit()

    if args.makeqa:
        print('Making QA Plots...')
        fastfit, metadata, coadd_type, _ = read_fastspecfit(fastfitfile)
        targetids = metadata['TARGETID'].data

        Spec.select(redrockfiles=redrockfile, targetids=targetids, use_quasarnet=False)
        data = Spec.read_and_unpack(CFit, fastphot=False, synthphot=True, remember_coadd=True)
        
        indx = np.arange(len(data))
        qaargs = [(CFit, EMFit, data[igal], fastfit[indx[igal]], metadata[indx[igal]],
                   coadd_type, False, qadir, None) for igal in np.arange(len(indx))]                

        if args.mp > 1:
            import multiprocessing
            with multiprocessing.Pool(args.mp) as P:
                P.map(_desiqa_one, qaargs)
        else:
            [desiqa_one(*_qaargs) for _qaargs in qaargs]
    else:
        Spec.select(redrockfiles=redrockfile, targetids=targetids, zmin=-0.1,
                    use_quasarnet=False, ntargets=len(zcat))
        data = Spec.read_and_unpack(CFit, fastphot=False, synthphot=True, remember_coadd=True)

        out, meta = Spec.init_output(CFit=CFit, EMFit=EMFit, fastphot=False)

        # Fit in parallel
        t0 = time.time()
        fitargs = [(iobj, data[iobj], out[iobj], meta[iobj], CFit, EMFit, False, False) # verbose and broadlinefit
                   for iobj in np.arange(Spec.ntargets)]
        if args.mp > 1:
            import multiprocessing
            with multiprocessing.Pool(args.mp) as P:
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

if __name__ == '__main__':
    main()
