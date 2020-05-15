#!/bin/bash

#SBATCH -n 3
#SBATCH -t 0-00:30
#SBATCH -o /home/gsvance/prismtools/output/prismMPI.%j.out
#SBATCH -e /home/gsvance/prismtools/output/prismMPI.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gsvance@asu.edu

# The python/2.7.9 module is necessary to run mpi4py
# I suspect the intel/2017x module will be needed to run PRISM itself
module purge
module load python/2.7.9
module load numpy/python-2x
#module load intel/2017x

TOOLSDIR="/home/gsvance/prismtools/"
SDFDIR="/home/gsvance/vconv_3d_asym_data/vconv_snsph_sdf_data/"
OUTDIR="/home/gsvance/vconv_3d_asym_data/vconv_prism_postprocess/"

cd ${TOOLSDIR}
mpiexec -n 3 python prism_mpi_wrapper.py ${SDFDIR} ${OUTDIR}

