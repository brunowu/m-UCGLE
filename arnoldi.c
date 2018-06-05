#include <stdio.h>
#include <mpi.h>
#include "Libs/mpi_lsa_com.h"

int main( int argc, char *argv[] ){
  MPI_Init( &argc, &argv );

  int arank, asize;

  int exit_type;

  MPI_Comm COMM_FATHER;

  MPI_Comm_size( MPI_COMM_WORLD, &asize );
  MPI_Comm_rank( MPI_COMM_WORLD, &arank );

  MPI_Comm_get_parent( &COMM_FATHER );

  if(arank == 0){
    printf("Info ]> The Comm world size of ERAM is %d \n", asize);
  }

  //check if any type to receive
  mpi_lsa_com_type_recv(&COMM_FATHER, &exit_type);

  //exit if receive the exit signal
  if(exit_type == 666){
    if(arank == 0){
      printf("Info ]> ERAM Receive quit information from Father and then exit\n");
    }
  }
  MPI_Finalize();

  return 0;

}
