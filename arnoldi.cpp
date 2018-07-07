#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"

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

  int length = 5;
  int i;

  double *data = new double [length];

  double set[5] = {1,2,3,4,5};

  for(i = 0; i < length; i++){
    data[i] = set[i];
  }

  mpi_lsa_com_array_send(&COMM_FATHER, &length, data);

  //check if any type to receive

  //check if any type to receive
  if(!mpi_lsa_com_type_recv(&COMM_FATHER, &exit_type)){
    if(arank == 0){
      printf("Info ]> ERAM Receive signal information from Father\n");
    }
  }

  //exit if receive the exit signal
  if(exit_type == 666){
    if(arank == 0){
      printf("Info ]> ERAM exit\n");
    }
  }

  free(data);

  MPI_Finalize();

  return 0;

}
