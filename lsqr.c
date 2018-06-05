#include <stdio.h>
#include <mpi.h>
#include "Libs/mpi_lsa_com.h"

int main( int argc, char *argv[] ){

  MPI_Init( &argc, &argv );

  int lrank, lsize;
  int exit_type;

  MPI_Comm COMM_FATHER;

  MPI_Comm_size( MPI_COMM_WORLD, &lsize );
  MPI_Comm_rank( MPI_COMM_WORLD, &lrank );

  MPI_Comm_get_parent( &COMM_FATHER );

  if(lrank == 0){
    printf("Info ]> The Comm world size of LS is %d \n", lsize);
  }

  //check if any type to receive
  mpi_lsa_com_type_recv(&COMM_FATHER, &exit_type);

  //exit if receive the exit signal
  if(exit_type == 666){
    if(lrank == 0){
      printf("Info ]> LS Receive quit information from Father and then exit\n");
    }
  }

  MPI_Finalize();

  return 0;

}
