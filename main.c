#include <stdio.h>
#include <mpi.h>
#include "Libs/mpi_lsa_com.h"

int main( int argc, char *argv[] ) {

  MPI_Init( &argc, &argv );

  int size, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if(rank == 0){
    printf("Info ]> The Comm world size of FATHER is %d \n", size);
  }

  int GMRES_SIZE = 2, ARNOLDI_SIZE = 2, LS_SIZE = 1;

  MPI_Comm COMM_GMRES, COMM_ARNOLDI, COMM_LS;
  MPI_Request gReq[GMRES_SIZE], aReq[ARNOLDI_SIZE], lReq[LS_SIZE];
  MPI_Status gStatus, aStatus, lStatus;

  MPI_Comm_spawn( "./gmres.exe", MPI_ARGV_NULL, GMRES_SIZE, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_GMRES, MPI_ERRCODES_IGNORE);
  MPI_Comm_spawn( "./arnoldi.exe", MPI_ARGV_NULL, ARNOLDI_SIZE, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_ARNOLDI, MPI_ERRCODES_IGNORE);
  MPI_Comm_spawn( "./lsqr.exe", MPI_ARGV_NULL, LS_SIZE, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_LS, MPI_ERRCODES_IGNORE);

  int exit_type;

  //receive exit type from GMRES Componet
  mpi_lsa_com_type_recv(&COMM_GMRES, &exit_type);

  int out_sended_type_a = 0, out_sended_type_l = 0;
  //send exit type to LS and ERAM Components
  mpi_lsa_com_type_send(&COMM_ARNOLDI, &exit_type, &out_sended_type_a);
  mpi_lsa_com_type_send(&COMM_LS, &exit_type, &out_sended_type_l);

  printf("Info ]> Father send exit type to ERAM and LS Component\n");
  printf("Debug ]> out_sended_type_a = %d / out_sended_type_l = %d \n", out_sended_type_a, out_sended_type_l);

  double *data;
  int length = 5;

  data = (double *)malloc(length*sizeof(double));
  double set[5] = {1,2,3,4,5};
  
  for(int i = 0; i < length; i++){
    data[i] = set[i];
  }

  mpi_lsa_com_array_send(&COMM_GMRES, &length, data);

  for(int i = 0; i < length; i++){
    printf("Debug ]>: Father send data[%d] = %f\n", i, data[i] );
  }

  free(data);

  MPI_Finalize();

  return 0;
}
