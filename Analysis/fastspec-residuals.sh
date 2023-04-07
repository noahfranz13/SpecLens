#!/bin/bash

source /global/cfs/cdirs/desi/software/desi_environment.sh 23.1
module swap desispec/0.57.0
#module load fastspecfit/2.1.1

codedir=$HOME
for package in fastspecfit; do
    echo Loading local check-out of $package
    export PATH=$codedir/$package/bin:$PATH
    export PYTHONPATH=$codedir/$package/py:$PYTHONPATH
done

dir=/global/cfs/cdirs/desi/users/nrfran/speclens
rrfile=$dir/redrock-source.fits
out=$dir/fastspec-source.fits
qa=$dir/QA-source

# run fastspecfit
fastspec $rrfile --mp $1 --outfile $out --specproddir custom 

# build QA plots
fastspecfit-qa $out -o $qa --redrockfiles $rrfile --mp $1
python3 displayQA.py --dir $qa
