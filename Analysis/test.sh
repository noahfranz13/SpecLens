#/bin/bash

source /global/cfs/cdirs/desi/software/desi_environment.sh 23.1
module load fastspecfit/2.1.1

#codedir=$HOME
#for package in fastspecfit; do
#    echo Loading local check-out of $package
#    export PATH=$codedir/$package/bin:$PATH
#    export PYTHONPATH=$codedir/$package/py:$PYTHONPATH
#done

python3 test.py