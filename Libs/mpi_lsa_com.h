#ifndef _MPI_LSA_COM_H_
#define _MPI_LSA_COM_H_

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

/*send type information functionality*/
int mpi_lsa_com_type_recv(MPI_Comm * com, int * type);
/*receive type information functionality*/
int mpi_lsa_com_type_send(MPI_Comm * com, int * type, int *count);

/*send array functionality*/
int mpi_lsa_com_array_send(MPI_Comm * com, int * size, double * data);

/*receive array functionality*/
int mpi_lsa_com_array_recv(MPI_Comm * com, int * size, double * data);

#endif
