#ifndef _MPI_LSA_COM_H_
#define _MPI_LSA_COM_H_

#include <stdio.h>
#include <mpi.h>

/*send type information functionality*/
int mpi_lsa_com_type_recv(MPI_Comm * com, int * type);
/*receive type information functionality*/
int mpi_lsa_com_type_send(MPI_Comm * com, int * type);

#endif
