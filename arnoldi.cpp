#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"
#include <complex>
#include <unistd.h>


//Block KrylovSchur METHODs to approximate the eigevalues

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziTpetraAdapter.hpp" //Anasazi interface to Tpetra
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include <Teuchos_CommandLineProcessor.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>

using Tpetra::CrsMatrix;
using Tpetra::Map;
using std::vector;


int main( int argc, char *argv[] ){
  MPI_Init( &argc, &argv );


  int arank, asize;

  int exit_type = 0;

  MPI_Comm COMM_FATHER;

  MPI_Comm_size( MPI_COMM_WORLD, &asize );
  MPI_Comm_rank( MPI_COMM_WORLD, &arank );

  MPI_Comm_get_parent( &COMM_FATHER );

  if(arank == 0){
    printf("Info ]> The Comm world size of ERAM is %d \n", asize);
  }

  int length = 5;
  int i, end = 0, p = 0;

  std::complex<double> *data = new std::complex<double> [length];

  double real[5] = {-0.830999, -0.774165, -0.463353, -0.444916, -0.287419};
  double imag[5] = {0.514128, 0.424414, 0.363388, 0.517988, 0.356902};

  while(!end){


    //check if any type to receive
    if(!mpi_lsa_com_type_recv(&COMM_FATHER, &exit_type)){
      if(arank == 0){
        printf("Info ]> ERAM Receive signal information from Father\n");
      }
      //exit if receive the exit signal
      if(exit_type == 666){
        if(arank == 0){
          printf("Info ]> ERAM exit\n");
        }
        end = 1;
        break;
      }
    }

    for(i = 0; i < length; i++){
      data[i].real(real[i]);
      data[i].imag(imag[i]);
    }

    usleep(100);
    mpi_lsa_com_cplx_array_send(&COMM_FATHER, &length, data);
    p++;
    printf("p = %d\n", p);
  }

  delete [] data;

  MPI_Comm_free(&COMM_FATHER);

  MPI_Finalize();

  return 0;

}
