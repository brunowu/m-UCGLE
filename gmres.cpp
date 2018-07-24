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
  int exit_type = 0;

  MPI_Comm COMM_FATHER;

  MPI_Comm_size( MPI_COMM_WORLD, &gsize );
  MPI_Comm_rank( MPI_COMM_WORLD, &grank );

  MPI_Comm_get_parent( &COMM_FATHER );

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("GMRES ]> Hello world from processor %s, rank %d out of %d processors\n", processor_name, grank, gsize);


/*
  if(grank == 0){
    printf("GMRES ]> The Comm world size of GMRES is %d \n", gsize);
  }
*/
  int length = EIGEN_ALL;
  int i;

  double *data_recv = new double [length];
  int end = 0;

  int count = 0;

  while(!end){
    if(!mpi_lsa_com_array_recv(&COMM_FATHER, &length, data_recv)){
      for(i = 0; i < length; i++){
        printf("GMRES ]>: GMRES rank = %d, data[%d] = %f\n",grank, i, data_recv[i] );
      }

      usleep(10000);
    } else{
      usleep(10000);
        //printf("GMRES NOT RECEVIECE ARRAY FROM FATHER\n");
    }

    count ++ ;

    if(count == 2000){
      end = 1;
      break;
    }

  }

  mpi_lsa_com_type_send(&COMM_FATHER, &type);

  if(grank == 0){
    printf("GMRES ]> GMRES send exit signal\n");
  }

  int gmres_final_exit;

  if(!mpi_lsa_com_type_recv(&COMM_FATHER, &gmres_final_exit)){
    if(gmres_final_exit == 777){
      printf("GMRES is allowed to exit by FATHER now\n");
    }else{
      usleep(1000000);
    }
  }else{
    usleep(10000000);
  }

  free(data_recv);

  MPI_Comm_free(&COMM_FATHER);

  printf("GMRES ]> Close of GMRES after waiting a little instant\n");
  MPI_Finalize();


  return 0;

}
