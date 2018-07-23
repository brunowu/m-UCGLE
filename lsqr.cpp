#include <stdio.h>
#include <cstdio>
#include <mpi.h>
#include "Libs/mpi_lsa_com.hpp"
#include <complex>
#include "Solvers/Utils/libs.hpp"
#include "Solvers/Utils/convhull.hpp"
#include "Solvers/Utils/ellipse.hpp"
#include "Solvers/LS/precond.hpp"
#include <cstring>
#include <iostream>
#include <complex>
#include <cmath>
#include <unistd.h>


#ifndef EIGEN_MIN
#define EIGEN_MIN 2
#endif

#ifndef EIGEN_ALL
#define EIGEN_ALL 10
#endif

#ifndef EIGEN_MAX
#define EIGEN_MAX 10
#endif

#ifndef MAX_PATH_LEN
#define MAX_PATH_LEN 4096
#endif

int main( int argc, char *argv[] ){

  MPI_Init( &argc, &argv );

  int lrank, lsize;
  int exit_type;

  MPI_Comm COMM_FATHER;

  MPI_Comm_size( MPI_COMM_WORLD, &lsize );
  MPI_Comm_rank( MPI_COMM_WORLD, &lrank );

  MPI_Comm_get_parent( &COMM_FATHER );

  if(lrank == 0){
    printf("Info ]> The Comm world size of LS is %d \n", lsize);
  }

  /* variables */
  int end, cumul, eigen_received, eigen_total, eigen_max;
  char load_path[MAX_PATH_LEN], export_path[MAX_PATH_LEN];
  int i, info, type = 0;

  int vector_size = 20;
  /*flags*/
  int flag, data_load = 0, data_export, continuous_export, data_load_any = 0;

  std::complex<double> *data, *eigen_cumul, *eigen_tri, *d, *c;

  double a_ell, c_ell, d_ell, d_reel;

  int data_size;
  int chsign;
  int mu1, mu2, mu, result_array_size;
  double *eta, *beta, *delta, alpha, scalar_tmp, *result_array;
  int ls_eigen_min, ls_eigen; // use default values

  sprintf(load_path,"./lsqr.bin");
  sprintf(export_path,"./lsqr.bin");

  /*options and flags we add later using Trilinos functions*/
  ls_eigen_min=EIGEN_MIN;
  ls_eigen=EIGEN_ALL;
  eigen_max=ls_eigen;

  eigen_tri = new std::complex<double> [vector_size];
  eigen_cumul = new std::complex<double> [vector_size];
  d = new std::complex<double> [vector_size + 1];
  c = new std::complex<double> [vector_size + 1];
  data = new std::complex<double> [vector_size];

  /* data that will be sended to GMRES for it's preconditionning step */
  eta = new double [eigen_max + 1];
  beta = new double [eigen_max + 1];
  delta = new double [eigen_max + 1];
  result_array = new double [eigen_max * 3 + 2];
  result_array_size = 2 + 3 * eigen_max;

  for(i = 0; i < (vector_size); i++){
    eigen_tri[i].real(0.0);
    eigen_cumul[i].real(0.0);

    eigen_tri[i].imag(0.0);
    eigen_cumul[i].imag(0.0);
  }

  for(i = 0;i < (vector_size) + 1; i++){
    d[i].real(0.0);
    c[i].real(0.0);

    d[i].imag(0.0);
    c[i].imag(0.0);
  }

  cumul = 0;
  eigen_received = 0;
  eigen_total = 0;
  end = 0;
  eigen_total = 0;
  ls_eigen = 0;

  int length = 5;

  while(!end){

    //exit type receiving and sending operation
    //check if any type to receive
    if(!mpi_lsa_com_type_recv(&COMM_FATHER, &exit_type)){
          printf("faieieieieie\n");
      printf("exit type of LS = %d\n", exit_type);
      //exit if receive the exit signal
      if(exit_type == 666){
        if(lrank == 0){
          printf("Info ]> LS Receive signal information from Father\n");
          printf("Info ]> LS exit\n");
        }
        end  = 1;
        break;
      }
    }
        /*in any case clear data array*/
    for(i = 0;i < eigen_max; i++){
      data[i].real(0.0);
      data[i].imag(0.0);
    }

    if(!mpi_lsa_com_cplx_array_recv(&COMM_FATHER, &data_size, data)){
      /* we received data or load it depending on the flags (for first step only*/

      for(i = 0; i < data_size; i++){
        printf("data[%d] = %f+%fi\n",i, data[i].real(),data[i].imag());
      }

      /* first we gonna remove some non-needed values */
      epurer(data,&data_size);

      for(i = 0;i < data_size; i++){
        printf("data_epure[%d] = %f+i%f\n",i,data[i].real(), data[i].imag() );
      }

      /* add them to the accumulated eigenvalues */
      /* if full renew full eigenvalues */
      if(eigen_total + data_size > vector_size) eigen_total = 0;
      /* select eigenvalues */
      for(i = 0; i < data_size; i++){
        eigen_cumul[eigen_total + i].real(data[i].real());
        eigen_cumul[eigen_total + i].imag(data[i].imag());
      }

      eigen_total += data_size;

      if(cumul < eigen_total) cumul = eigen_total;

      for(i = 0;i < cumul; i++){

        eigen_tri[i].real(eigen_cumul[i].real());
        eigen_tri[i].imag(eigen_cumul[i].imag());

      }
      
      eigen_received += data_size;
      /* if we didn't received enough eigenvalues */
      if(eigen_received < ls_eigen_min && data){
          printf("eigen_received = %d\n",eigen_received );
          continue;
      }
      else{
        eigen_received = 0;
        tri(eigen_tri, cumul, &chsign);

        mu1=0;

        mu2=0;

        printf("FUCKKKKKKKK  YYYYOUUUUUUUUU\n");


        /* convex hull computation */
        if(chsign > 0){
        //keepPositif(eigen_tri,&cumul);
          convhull(eigen_tri, c, d, chsign, &mu1, 0, 0);
          printf("@} LSQR convhul negatif chsigne %d cumul %d mu1 %d\n",chsign,cumul,mu1);
        }
        if(chsign < cumul){
          convhull(eigen_tri, c, d, cumul-chsign, &mu2, chsign, mu1);
          printf("@} LSQR convhul positif chsigne %d cumul %d mu1 %d mu2 %d\n",chsign,cumul,mu1,mu2);
        }
        mu = mu1 + mu2;

        /* Ellipse computation */
        ellipse(c,  d, mu+1, mu1, &c_ell, &a_ell, &d_ell, &d_reel, &info);

        if(fabs(d_ell)<epsilon()) d_ell = 1.;

        if(fabs(a_ell)<epsilon()){
          ls_eigen=0;
        }
        else{
          LSPrecond(a_ell, d_ell,c_ell,eta, &alpha, beta,
                delta, c, d,&mu, &ls_eigen, &ls_eigen_min, &eigen_max); 
        }

        for(i = 0; i < eigen_max; i++){
          //printf("eta[%d] = %f, alpha = %f, beta[%d] = %f, deta[%d] =  %f\n", i, eta[i], alpha, i, beta[i],i, delta[i]);
        }

      }
    }

    if(ls_eigen > 1){
      /* place the computed results inside the array */
      scalar_tmp = ls_eigen;
      std::memcpy(&result_array[0],&scalar_tmp,1*sizeof(double));
      std::memcpy(&result_array[1],&alpha,1*sizeof(double));
      std::memcpy(&result_array[2],eta,ls_eigen*sizeof(double));
      std::memcpy(&result_array[2 + ls_eigen],beta,ls_eigen*sizeof(double));
      std::memcpy(&result_array[2 + 2 * ls_eigen],delta,ls_eigen*sizeof(double));

      /*
        if(continuous_export){
          to do: data export
        }
        */
      /* and send it */
      result_array_size=2+3*ls_eigen;

      MPI_Request rq[1];
      mpi_lsa_com_array_send(&COMM_FATHER, &result_array_size,result_array, rq);
      printf("LS ]> LS send array\n");
    }

    printf("didudiu\n");

    if(ls_eigen>1){
      ls_eigen=0;
    }
  }

  /*
  if(data_export){
    to do: data export
  }
  */

  delete [] eigen_tri;
  delete [] eigen_cumul;
  delete [] d;
  delete [] c;
  delete [] data;

  delete [] eta;
  delete [] beta;
  delete [] delta;
  delete [] result_array;

  MPI_Comm_free(&COMM_FATHER);
  MPI_Finalize();

  return 0;

}
