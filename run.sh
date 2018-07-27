#!/bin/bash

#mpirun -np 1 ./test.exe -nb_gmres 2 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe ./gmres2.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe

#mpirun -np 1 ./test.exe -nb_gmres 1 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe --filename=\"utm300.mtx\"

mpirun -np 1 ./test.exe -nb_gmres 1 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe --filename=\"utm300.mtx\" --eps-quiet --eps-nodebug --eps-exsitu --eps-sort=\"LM\" \ --eps-nonherm --eps-nev=10 --eps-blockSize=1 --eps-tol=0.1 --eps-no-print-matrix --eps-all-print --eps-numBlocks=100 --eps-maxRestarts=100 --ksp-quiet --ksp-nodebug --ksp-frequency=-1 --ksp-tol=1e-05 --ksp-num-rhs=2 --ksp-block-size=100 --ksp-no-precond --ksp-num-blocks=50 --ksp-fixed-tol --ksp-no-print-matrix --ksp-all-print --ksp-no-dump-data --lsp-use