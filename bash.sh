#!/bin/bash

cd Libs
mpicc -c mpi_lsa_com.c
cd ..
mpicc main.c -o test.exe Libs/mpi_lsa_com.o
mpicc arnoldi.c -o arnoldi.exe Libs/mpi_lsa_com.o
mpicc arnoldi.c -o arnoldi2.exe Libs/mpi_lsa_com.o
mpicc gmres.c -o gmres.exe Libs/mpi_lsa_com.o
mpicc gmres.c -o gmres2.exe Libs/mpi_lsa_com.o
mpicc lsqr.c -o lsqr.exe Libs/mpi_lsa_com.o

mpirun -n 1 ./test.exe

