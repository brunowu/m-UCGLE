#include <stdio.h>
#include <mpi.h>

int main( int argc, char *argv[] ){
  MPI_Init( &argc, &argv );

  int grank, gsize;

  MPI_Comm_size( MPI_COMM_WORLD, &gsize );
  MPI_Comm_rank( MPI_COMM_WORLD, &grank );

  if(grank == 0){
    printf("Info ]> The Comm world size of LS is %d \n", gsize);
  }

  MPI_Finalize();

  return 0;

}
