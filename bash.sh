#!/bin/bash

cd Libs
mpicc -c mpi_lsa_com.c
mpicc -c args_parser.c

cd ..
mpicc main.c -o test.exe Libs/mpi_lsa_com.o Libs/args_parser.o
mpicc arnoldi.c -o arnoldi.exe Libs/mpi_lsa_com.o Libs/args_parser.o
mpicc arnoldi.c -o arnoldi2.exe Libs/mpi_lsa_com.o Libs/args_parser.o
mpicc gmres.c -o gmres.exe Libs/mpi_lsa_com.o Libs/args_parser.o
mpicc gmres.c -o gmres2.exe Libs/mpi_lsa_com.o Libs/args_parser.o
mpicc lsqr.c -o lsqr.exe Libs/mpi_lsa_com.o Libs/args_parser.o

mpirun -np 1 ./test.exe -nb_gmres 2 -proc_gmres 1 -nb_arnoldi 2 -proc_arnoldi 1 -gmres_exec ./gmres.exe ./gmres2.exe -arnoldi_exec ./arnoldi.exe ./arnoldi2.exe -lsqr_exec ./lsqr.exe
