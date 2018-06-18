#!/bin/bash
#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "TEST 1"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=02:30:00

#SBATCH -n 8

mpirun -np 8 ./blockGMRES --verbose --tol=1e-06 --block-size=2 --num-blocks=200 --frequency=10 --num-rhs=2 #--filename="mhd1280b.cua"


