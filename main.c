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

  double *data;
  int length = 5;
  int i;

  double *data_recv;

  data = (double *)malloc(length*sizeof(double));
  data_recv = (double *)malloc(length*sizeof(double));

  //receive data from ERAM
  if(!mpi_lsa_com_array_recv(&COMM_ARNOLDI, &length, data_recv)){
    printf("Info ]> Array receive from ERAM Component\n");
  }

  //send this array to LS
  mpi_lsa_com_array_send(&COMM_LS, &length, data_recv);
  printf("Info ]> Father send array to LS\n");
  //receive new array from LS
  if(!mpi_lsa_com_array_recv(&COMM_LS, &length, data)){
    printf("Info ]> Father has Array received from LS Component\n" );
    for(i = 0; i < length; i++){
      printf("Debug ]>: Father send data[%d] = %f to GMRES\n", i, data[i] );
    }
    //send new array to GMRES
    mpi_lsa_com_array_send(&COMM_GMRES, &length, data);
  }

  int exit_type;

  //receive exit type from GMRES Componet
  if(!mpi_lsa_com_type_recv(&COMM_GMRES, &exit_type)){
    printf("Info ]> Father Receive exit type from GMRES Component\n" );
  }

  int out_sended_type_a = 0, out_sended_type_l = 0;
  //send exit type to LS and ERAM Components
  mpi_lsa_com_type_send(&COMM_ARNOLDI, &exit_type, &out_sended_type_a);
  mpi_lsa_com_type_send(&COMM_LS, &exit_type, &out_sended_type_l);

  printf("Info ]> Father send exit type to ERAM and LS Component\n");
  printf("Debug ]> out_sended_type_a = %d / out_sended_type_l = %d \n", out_sended_type_a, out_sended_type_l);

  free(data);
  free(data_recv);

  MPI_Finalize();

  return 0;
}
