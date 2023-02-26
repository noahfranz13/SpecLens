#!/bin/bash

source /global/cfs/cdirs/desi/software/desi_environment.sh 23.1
module load fastspecfit/2.1.1

dir=/global/cfs/cdirs/desi/users/nrfran/speclens
file=$dir/fastspec-input.fits

# crossmatch lenses
cd FindLenses
python3 crossmatch_known_lenses.py --specprod iron --outdir $dir

# separate found lenses
cd ../Analysis
python3 runSpecLens.py -i $file --makeqa 
