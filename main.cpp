#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"
#include "Libs/args_parser.hpp"
#include "utils/logo.h"
#include <complex>
#include <unistd.h>

#ifndef EIGEN_ALL
#define EIGEN_ALL 1000
#endif


int main( int argc, char *argv[] ) {

  int gmres_nb, arnoldi_nb, gmres_proc, arnoldi_proc, ls_proc = 1;
  char **gmres_cmds;
  char **arnoldi_cmds;
  char *lsqr_cmd;

  MPI_Init( &argc, &argv );

  int size, rank;
  int i, j, end = 0;

  logo();

  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //Main world of MUCGLE should be 1
  if(size != 1){
    printf("Error ]> The Comm world size of manager engine should be 1\n");
    return 1;
  }

  border_print2();
  center_print("Start Initialization", 79);
  border_print2();
  //Parser the arguments to initialize MUCLGE
  argsParser(argc, argv, &gmres_nb, &arnoldi_nb, &gmres_proc, &arnoldi_proc);

  //Parser the executables string name for spawning
  gmres_cmds = argsParserGmresExec(argc, argv, gmres_nb);
  arnoldi_cmds = argsParserArnoldiExec(argc, argv, arnoldi_nb);
  lsqr_cmd = argsParserLsqrExec(argc, argv);

  MPI_Comm COMM_GMRES[2], COMM_ARNOLDI[2], COMM_LS;
  MPI_Request gReq[gmres_proc], aReq[arnoldi_proc], lReq[ls_proc];
  MPI_Status gStatus, aStatus, lStatus;

  border_print2();
  center_print("Start Spawning Executables", 79);
  border_print2();

  for(i = 0; i < gmres_nb; i++){
    MPI_Comm_spawn( gmres_cmds[i], MPI_ARGV_NULL, gmres_proc, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_GMRES[i], MPI_ERRCODES_IGNORE);
  }

  for( j = 0; j < arnoldi_nb; j++){
    MPI_Comm_spawn( arnoldi_cmds[j], MPI_ARGV_NULL, arnoldi_proc, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_ARNOLDI[j], MPI_ERRCODES_IGNORE);
  }

  MPI_Comm_spawn( lsqr_cmd, MPI_ARGV_NULL, ls_proc, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &COMM_LS, MPI_ERRCODES_IGNORE);

  border_print2();
  center_print("Start Resolving Linear Systems by MUCGLE", 79);
  border_print2();

  int length = 5;

  int exit_type[gmres_nb];
  int exit_recv[gmres_nb];
  int exit_signal;
  int exit = 0;
  int exit_signal_gmres = 777;

  int out_sended_type_a[2] = {0,0};
  int out_sended_type_g[2] = {0,0};
  int out_sended_type_l;

  MPI_Request rq[gmres_nb];
  MPI_Status status[gmres_nb];
  int flg;

  std::complex<double> *data_recv = new std::complex<double> [length];
  double *data = new double [EIGEN_ALL];

  for(i = 0; i < gmres_nb; i++){
    exit_type[i] = 0;
  }
  int tmp;

  while(!end){

    
    //receive exit type from GMRES Componet
    for(i = 0; i < gmres_nb; i++){
      if(!mpi_lsa_com_type_recv(&COMM_GMRES[i], &tmp)){
        exit_type[i] = tmp;
        printf("tmp = %d\n", tmp);
        printf("Main ]> Father Receive exit type from GMRES %d Component = %d\n", i, exit_type[i]);
        if (exit_type[i] == 666){
          exit = exit + 1;
        }
      }
    }

    for(i = 0; i < gmres_nb; i++){
      printf("exit_type[%d] = %d\n", i, exit_type[i]);
    }

    if(exit == gmres_nb){
      exit_signal = 666;
      //send exit type to LS and ERAM Components
      printf("Main ]> Father send exit type to ERAM and LS Component\n");

      for(j = 0; j < arnoldi_nb; j++){
        mpi_lsa_com_type_send(&COMM_ARNOLDI[j], &exit_signal);
      }
      mpi_lsa_com_type_send(&COMM_LS, &exit_signal);
      end = 1;
      break;
    }
    
    //receive data from ERAM
    for(j = 0; j < arnoldi_nb; j++){
      if(!mpi_lsa_com_cplx_array_recv(&COMM_ARNOLDI[j], &length, data_recv)){
        printf("Main ]> Array receive from ERAM %d Component\n", j);
      }
    }

    //send this array to LS
    mpi_lsa_com_cplx_array_send(&COMM_LS, &length, data_recv);
    printf("Main ]> Father send array to LS\n");
    //receive new array from LS
    if(!mpi_lsa_com_array_recv(&COMM_LS, &length, data)){
      printf("Main ]> Father has Array received from LS Component\n" );
      for(i = 0; i < length; i++){
        //printf("Debug ]>: Father send data[%d] = %f to GMRES\n", i, data[i] );
      }
      //send new array to multiples GMRES
      for(i = 0; i < gmres_nb; i++){
        mpi_lsa_com_array_send(&COMM_GMRES[i], &length, data, rq);
        printf("Main ]> Father has sent preconditioned parameters to GMRES Component\n" );
      }
    }
  }

  for(i = 0; i < gmres_nb; i++){
    mpi_lsa_com_type_send(&COMM_ARNOLDI[i], &exit_signal_gmres);
    printf("Main ]> Father has final exit signal to GMRES Component\n" );
  }
  
  border_print2();
  center_print("Remove Application", 79);
  border_print2();

  delete [] data;
  delete [] data_recv;

  for(i = 0; i < gmres_nb; i++){
      delete [] gmres_cmds[i];
  }

  for(i = 0; i < arnoldi_nb; i++){
      delete [] arnoldi_cmds[i];
  }

  delete [] lsqr_cmd;

  for(j = 0; j < gmres_nb; j++){
    MPI_Comm_free(&COMM_GMRES[j]);
  }

  for(j = 0; j < arnoldi_nb; j++){
    MPI_Comm_free(&COMM_ARNOLDI[j]);
  }

  MPI_Comm_free(&COMM_LS);

  MPI_Finalize();

  return 0;
}
