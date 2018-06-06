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

  double *data_recv;
  int length = 5;
  int i;

  data_recv = (double *)malloc(length*sizeof(double));

  if(!mpi_lsa_com_array_recv(&COMM_FATHER, &length, data_recv)){
    printf("done\n");
    for(i = 0; i < length; i++){
      printf("Debug ]>: GMRES rank = %d, data[%d] = %f\n",grank, i, data_recv[i] );
    }
  }

  int out_sended_type = 0;
  mpi_lsa_com_type_send(&COMM_FATHER, &type, &out_sended_type);
  if(grank == 0){
    printf("Debug ]> out_sended_type = %d\n", out_sended_type);
  }

  free(data_recv);

  MPI_Finalize();

  return 0;

}
