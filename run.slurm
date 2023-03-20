#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --job-name=speclens
#SBATCH --account=desi
#SBATCH --mail-user=nr25fran@siena.edu
#SBATCH --mail-type=ALL
#SBATCH -o /global/cfs/cdirs/desi/users/nrfran/speclens/speclens.log
#SBATCH -c 128

#srun -c 128 ./run.sh 10
./run.sh 128
