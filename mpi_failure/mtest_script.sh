#!/bin/bash

#SBATCH -n 5
#SBATCH -t 0-00:03
#SBATCH -o /home/gsvance/prismtools/mt.%j.out
#SBATCH -e /home/gsvance/prismtools/mt.%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=gsvance@asu.edu

module purge
module load python/2.7.9

cd /home/gsvance/prismtools/
mpiexec -n 5 python mtest.py

