#!/bin/sh

mpicc args_parser.c 
mpirun -np 1 ./a.out -nb_gmres 1 -proc_gmres 2 -nb_arnoldi 3 -proc_arnoldi 4 -gmres_exec ./gmres1 -arnoldi_exec ./arnoldi1 ./arnoldi2 ./arnoldi3

