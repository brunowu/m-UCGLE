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

  double *data_recv;
  int length = 5;
  int i;
  double *data_send;

  data_recv = (double *)malloc(length*sizeof(double));
  data_send = (double *)malloc(length*sizeof(double));

  if(!mpi_lsa_com_array_recv(&COMM_FATHER, &length, data_recv)){
    for(i = 0; i < length; i++){
      data_send[i] = data_recv[i]*data_recv[i]+data_recv[i];
      printf("Debug ]>: LS rank = %d, data_recv[%d] = %f, data_send[%d] = %f\n",lrank, i, data_recv[i], i, data_send[i]);
    }
    mpi_lsa_com_array_send(&COMM_FATHER, &length, data_send);
  }

  //check if any type to receive
  if(!mpi_lsa_com_type_recv(&COMM_FATHER, &exit_type)){
    printf("Info ]> LS Receive signal information from Father\n");
  }

  //exit if receive the exit signal
  if(exit_type == 666){
    if(lrank == 0){
      printf("Info ]> LS exit\n");
    }
  }

  free(data_recv);
  free(data_send);

  MPI_Comm_free(&COMM_FATHER);

  MPI_Finalize();

  return 0;

}
