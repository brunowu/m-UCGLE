#!/bin/bash

#SBATCH --comment "Hello ROMEO!"
#SBATCH -J "TEST 1"

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=02:30:00

#SBATCH -n 8
#SBATCH -N 1

#mpirun -np 1 ./test.exe -nb_gmres 2 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe ./gmres2.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe

#mpirun -np 1 ./test.exe -nb_gmres 1 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe --filename=\"utm300.mtx\"

#module load gcc/7.2.0
#module load openmpi/2.1.1-gcc7.2.0

#mpiexec -n 1 ./test.exe -nb_gmres 1 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe --filename=\"mhd1280a.mtx\" --eps-quiet --eps-nodebug --eps-exsitu --eps-sort=\"LM\" --eps-nonherm --eps-nev=10 --eps-blockSize=1 --eps-tol=0.1 --eps-no-print-matrix --eps-all-print --eps-numBlocks=20 --eps-maxRestarts=50 --ksp-filename=\"mhd1280a.mtx\" --ksp-nodebug --ksp-frequency=10 --ksp-tol=1e-05 --ksp-num-rhs=2 --ksp-block-size=2 --ksp-no-precond --ksp-num-blocks=50 --ksp-fixed-tol --ksp-no-print-matrix --ksp-all-print --ksp-no-dump-data --ksp-lsp-degree=5 --ksp-lsp-latency=1 --ksp-use-lsp




mpirun -np 1 ./test.exe -nb_gmres 1 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe --filename=\"mhd1280a.mtx\" --eps-quiet --eps-nodebug --eps-exsitu --eps-sort=\"LM\" --eps-nonherm --eps-nev=10 --eps-blockSize=1 --eps-tol=0.1 --eps-no-print-matrix --eps-all-print --eps-numBlocks=20 --eps-maxRestarts=50 --ksp-filename=\"mhd1280a.mtx\" --ksp-nodebug --ksp-frequency=10 --ksp-tol=1e-05 --ksp-num-rhs=2 --ksp-block-size=2 --ksp-no-precond --ksp-num-blocks=50 --ksp-fixed-tol --ksp-no-print-matrix --ksp-all-print --ksp-no-dump-data --ksp-lsp-degree=5 --ksp-lsp-latency=1 --ksp-use-lsp
