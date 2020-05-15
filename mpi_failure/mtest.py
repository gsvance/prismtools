### !!!! Needs module load python/2.7.9 in order to run
### This is where mpi4py imports from!!!

from mpi4py import MPI
import sys


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#print rank

print rank, sys.argv

