#!/bin/bash

mpicc master.c -o master_program
mpicc worker.c -o worker_program

mpirun -n 1 master_program


