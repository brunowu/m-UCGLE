#!/bin/bash

cd Libs

cd ..
mpic++ main.cpp -o test.exe
mpic++ arnoldi.cpp -o arnoldi.exe
mpic++ arnoldi.cpp -o arnoldi2.exe 
mpic++ gmres.cpp -o gmres.exe
mpic++ gmres.cpp -o gmres2.exe
mpic++ lsqr.cpp -o lsqr.exe

mpirun -np 1 ./test.exe -nb_gmres 2 -proc_gmres 1 -nb_arnoldi 2 -proc_arnoldi 1 -gmres_exec ./gmres.exe ./gmres2.exe -arnoldi_exec ./arnoldi.exe ./arnoldi2.exe -lsqr_exec ./lsqr.exe
