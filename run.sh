#!/bin/bash

source /global/cfs/cdirs/desi/software/desi_environment.sh 23.1
module swap desispec/0.57.0
module load fastspecfit/2.1.1

dir=/global/cfs/cdirs/desi/users/nrfran/speclens
file=$dir/fastspec-input.fits
qa=$dir/QA

# crossmatch lenses
cd FindLenses
python3 crossmatch_known_lenses.py --specprod iron --outdir $dir

# separate found lenses
cd ../Analysis
python3 runSpecLens.py -i $file --makeqa --overwrite --mp $1 

# generate html page of QA plots
python3 displayQA.py --dir $qa