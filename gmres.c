#include <stdio.h>
#include <mpi.h>
#include "Libs/mpi_lsa_com.h"

int main( int argc, char *argv[] ){
  MPI_Init( &argc, &argv );

  int grank, gsize;
  int type = 666;

  MPI_Comm COMM_FATHER;

  MPI_Comm_size( MPI_COMM_WORLD, &gsize );
  MPI_Comm_rank( MPI_COMM_WORLD, &grank );

  MPI_Comm_get_parent( &COMM_FATHER );

  if(grank == 0){
    printf("Info ]> The Comm world size of GMRES is %d \n", gsize);
  }

  mpi_lsa_com_type_send(&COMM_FATHER, &type);

  MPI_Finalize();

  return 0;

}
