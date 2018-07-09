#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"
#include <complex>

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


  std::complex<double> *data = new std::complex<double> [length];

  double real[5] = {-0.830999, -0.774165, -0.463353, -0.444916, -0.287419};
  double imag[5] = {0.514128, 0.424414, 0.363388, 0.517988, 0.356902};

  for(i = 0; i < length; i++){
    data[i].real(real[i]);
    data[i].imag(imag[i]);
  }

  mpi_lsa_com_cplx_array_send(&COMM_FATHER, &length, data);

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

  delete [] data;

  MPI_Comm_free(&COMM_FATHER);

  MPI_Finalize();

  return 0;

}
