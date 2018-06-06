#include <stdio.h>
#include <mpi.h>
#include "Libs/mpi_lsa_com.h"

int main( int argc, char *argv[] ) {

  MPI_Init( &argc, &argv );

  int size, rank;
  int i, j;

  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  int gmres_nb = 2; //gmres number = spawn

  int arnoldi_nb = 2;

  char gmres_cmds[][20] = {"./gmres.exe", "./gmres2.exe"};

  char arnoldi_cmds[][20] = {"./arnoldi.exe", "./arnoldi2.exe"};

//  gmres_cmd = (char *) malloc (gmres_nb * sizeof (char));

  if(rank == 0){
    printf("Info ]> The Comm world size of FATHER is %d \n", size);
  }

  int GMRES_SIZE = 2, ARNOLDI_SIZE = 2, LS_SIZE = 1;

  MPI_Comm COMM_GMRES[2], COMM_ARNOLDI[2], COMM_LS;
  MPI_Request gReq[GMRES_SIZE], aReq[ARNOLDI_SIZE], lReq[LS_SIZE];
  MPI_Status gStatus, aStatus, lStatus;

  for(i = 0; i < gmres_nb; i++){
    MPI_Comm_spawn( gmres_cmds[i], MPI_ARGV_NULL, GMRES_SIZE, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_GMRES[i], MPI_ERRCODES_IGNORE);
  }
//  MPI_Comm_spawn( "./gmres2.exe", MPI_ARGV_NULL, GMRES_SIZE, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_GMRES2, MPI_ERRCODES_IGNORE);

  for( j = 0; j < arnoldi_nb; j++){
    MPI_Comm_spawn( arnoldi_cmds[j], MPI_ARGV_NULL, ARNOLDI_SIZE, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_ARNOLDI[j], MPI_ERRCODES_IGNORE);
  }

  MPI_Comm_spawn( "./lsqr.exe", MPI_ARGV_NULL, LS_SIZE, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_LS, MPI_ERRCODES_IGNORE);

  double *data;
  int length = 5;

  double *data_recv;

  data = (double *)malloc(length*sizeof(double));
  data_recv = (double *)malloc(length*sizeof(double));

  //receive data from ERAM
  for(j = 0; j < arnoldi_nb; j++){
    if(!mpi_lsa_com_array_recv(&COMM_ARNOLDI[j], &length, data_recv)){
      printf("Info ]> Array receive from ERAM %d Component\n", j);
    }
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
    //send new array to multiples GMRES
    for(i = 0; i < gmres_nb; i++){
      mpi_lsa_com_array_send(&COMM_GMRES[i], &length, data);
    }
  }

  int exit_type[gmres_nb];
  int exit_signal;
  int exit = 1;

  //receive exit type from GMRES Componet
  for(i = 0; i < gmres_nb; i++){
    if(!mpi_lsa_com_type_recv(&COMM_GMRES[i], &exit_type[i])){
      printf("Info ]> Father Receive exit type from GMRES %d Component\n", i );
    }
    exit &= (exit_type[i] == 666);
  }

  int out_sended_type_a[2] = {0,0};
  int out_sended_type_l;

  if(exit){
    exit_signal = 666;
    //send exit type to LS and ERAM Components
    for(j = 0; j < arnoldi_nb; j++){
      mpi_lsa_com_type_send(&COMM_ARNOLDI[j], &exit_signal, &out_sended_type_a[j]);
    }
    mpi_lsa_com_type_send(&COMM_LS, &exit_signal, &out_sended_type_l);

    printf("Info ]> Father send exit type to ERAM and LS Component\n");
  }

  free(data);
  free(data_recv);

  MPI_Finalize();

  return 0;
}
