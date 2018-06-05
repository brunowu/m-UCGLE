#include <stdio.h>
#include <mpi.h>
#include "../mpi_lsa_com.h"

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

  int flag = 0, count;
  MPI_Status status;
  MPI_Request request;

  while(!flag){
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,COMM_FATHER, &flag, &status);
  }

  if(flag){
    MPI_Get_count(&status,MPI_INT,&count);

    printf("Debug ]> Rank %d on ERAM: size = %d \n",arank, count );

    if(count == 1){
      MPI_Recv(&exit_type, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, COMM_FATHER, &status);
    }

    printf("Debug ]> Rank %d on ERAM: receive exit type = %d \n", arank, exit_type);
  }

  if(exit_type == 666){
    if(arank == 0){
      printf("Info ]> ERAM Receive quit information and then exit\n");
    }
  }
  MPI_Finalize();

  return 0;

}
