#!/bin/bash

#mpirun -np 1 ./test.exe -nb_gmres 2 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe ./gmres2.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe

mpirun -np 1 ./test.exe -nb_gmres 1 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe --filename=\"utm300.mtx\"

