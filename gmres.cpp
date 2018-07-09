#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"
#include <unistd.h>

#ifndef EIGEN_ALL
#define EIGEN_ALL 1000
#endif

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

  int length = EIGEN_ALL;
  int i;

  double *data_recv = new double [length];

  if(!mpi_lsa_com_array_recv(&COMM_FATHER, &length, data_recv)){
    for(i = 0; i < length; i++){
      //printf("Debug ]>: GMRES rank = %d, data[%d] = %f\n",grank, i, data_recv[i] );
    }
  }

  usleep(1000000);

  int out_sended_type = 0;
  mpi_lsa_com_type_send(&COMM_FATHER, &type, &out_sended_type);
  if(grank == 0){
    printf("Debug ]> GMRES out_sended_type = %d\n", out_sended_type);
  }

  free(data_recv);

  MPI_Comm_free(&COMM_FATHER);

  MPI_Finalize();

  return 0;

}
