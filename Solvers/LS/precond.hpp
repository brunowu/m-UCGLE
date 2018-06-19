#include "../Utils/libs.hpp"
#include <complex>
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Version.hpp"

int LSPrecond(double a_ell, double d_ell, double c_ell, std::complex<double> * eta, std::complex<double> * alpha, std::complex<double> * beta,
              std::complex<double> * delta, std::complex<double> * c, std::complex<double> * d, int * mu, int * nb_eigen, int * min_eigen,
              int * nb_eigen_all){

  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  // Creating an instance of the LAPACK class for double-precision routines looks like:
  Teuchos::LAPACK<int, double> lapack;

  //////////////////////////////
  std::complex<double> * fact_tmp;
  std::complex<double> * res;

  int i, j, k, nu;

  /* allocate work array gamma and mm_tmp, init it to 0 */
  std::complex<double> ** gamma  = new std::complex<double> * [(*nb_eigen_all)+3];

  for(i = 0; i < (*nb_eigen_all)+3; i++){
    gamma[i] = new std::complex<double> [(*nb_eigen_all)+3];
    for(j = 0; j < (*nb_eigen_all)+3; j++){
      gamma[i][j](0.0,0.0);
    }
  }

  std::complex<double> ** mm_tmp = new std::complex<double> * [(*nb_eigen_all)+1];

  for(i = 0; i < (*nb_eigen_all)+1; i++){
    mm_tmp[i] = new std::complex<double> [(*nb_eigen_all)+1];
    for(j = 0; j < (*nb_eigen_all)+1; j++){
      mm_tmp[i][j](0.0, 0.0);
    }
  }

  for(i = 0; i < *nb_eigen_all; i++){
    beta[i](0.0, 0.0);
    delta[i](0.0, 0.0);
  }

  /* begin computations */
  beta[i] = a_ell / 2.0;
  *nb_eigen = *nb_eigen_all;

  for(i = 0; i < *nb_eigen_all; i++){
    delta[i + 1] = (d_ell) / (4 * beta[i]);
    beta[i + 1] = (a_ell) - delta[i + 1];
    if(std::norm(beta[i + 1]) < epsilon()){
      *nb_eigen = i;
      if(*nb_eigen < * min_eigen){
        return 0;
      }
    }
  }

  *alpha = c_ell;

  /* computation of gamma */
  for(nu = 0; nu < *mu; nu++){
    for(i = 0; i < * nb_eigen_all + 1; i++){
      for(j =0; j < * nb_eigen_all + 1; j++){
        gamma[i][j] (0.0,0.0);
      }
    }
    gamma[1][1](1.0,0.0);

    for(j = 1; j < *nb_eigen_all; j++){
      for(i = 1; i < j + 1; i++){
        gamma[i][j + 1](d[nu].real()/2.0 * gamma[i + 1][j].real() + gamma[i - 1][j].real()
                        -(d[nu].imag()/2.0 + gamma[i + 1][j].imag()) + gamma[i - 1][j].imag()
                        +((c[nu]) - *alpha).real() * gamma[i][j].real() - (c[nu]).imag() * gamma[i][j].imag()
                        -delta[j - 1].real() * (gamma[i][j - 1].real()) / (beta[j - 1]).real(),

                        gamma[i][j + 1].imag());

      gamma[i][j + 1](gamma[i][j + 1].real(),

                      d[nu].real() / 2.0 * (gamma[i + 1][j].imag() + gamma[i - 1][j].imag())
                      + d[nu].imag() / 2.0 * (gamma[i + 1][j]).real() + gamma[i - 1][j].real()
                      + (c[nu] - *alpha).real() * (gamma[i][j]).imag()
                      + (c[nu].imag() * gamma[i][j].real())
                      - (delta[j - 1].real() + gamma[i][j - 1].imag() / beta[j - 1]));
      }
      gamma[0][j + 1] = gamma[2][j + 1];
    }

    /* computation of MM */
    for(j = 0; j <= *nb_eigen_all; j++){
      for(i = 0; i <= j; i++){
        mm_tmp[i][j] = mm_tmp[i][j] + 4.*(((gamma[1][j + 1].real())*(gamma[1][i + 1].real())
							    +(gamma[1][j + 1].imag())*(gamma[1][i + 1]).imag()));
        if(i > 1){
  					for(k = 1; k < i; k++){
  						mm_tmp[i][j] = mm_tmp[i][j] + 2.*(((gamma[k + 1][j + 1].real())*(gamma[k + 1][i + 1].real())
  						+(gamma[k + 1][j + 1].imag())*(gamma[k + 1][i + 1].imag())));
  					}
  			 }
      }
    }
  }

  Teuchos::SerialDenseMatrix<int, double> MM((*nb_eigen_all) + 1,(*nb_eigen_all) + 1);

  /* Filling of the lower triangular part */
  for(j = 0; j <= *nb_eigen_all; j++){
    for(i = 0; i <= j; i++){
      MM(i, j) = mm_tmp[i][j].real();
      MM(j, i) = mm_tmp[i][j].real();
    }
  }

  /*proceed to Cholesky factorization by the lapack interface provided by Teuchos*/
  int info;
  char UPLO = 'L';
  lapack.POTRF(UPLO, 3, MM.values(), MM.stride(), &info);

  Teuchos::SerialDenseMatrix<int, double> fact((*nb_eigen_all) + 1,(*nb_eigen_all) + 1);

  /*get Cholesky factorized matrix fact*/

  for(j = 0; j <= *nb_eigen_all; j++){
    for(i = 0; i <= j; i++){
      fact(i, j) = MM[i][j];
    }
  }

  /*Create the matrix operator F that will be used in the QR factorization*/
  Teuchos::SerialDenseMatrix<int, double> F((*nb_eigen) + 1,(*nb_eigen) + 1);
  F(0, 0) = (*alpha) * (fact.values())[0] + beta[0] * (fact.values())[1];
  F(1, 0) = beta[0] * (fact.values())[1 + (( *nb_eigen_all) + 1)];

  for(j = 1; j < *nb_eigen; j++){
    for(i = 0; i < j; i++){
      F(i, j) = delta[j-1]*(fact.values())[(j-1) + i*((*nb_eigen_all) + 1)]
			       + (*alpha) * (fact.values())[j + i * ((*nb_eigen_all) + 1)]
			       + (beta[j])* (fact.values())[(j + 1) + i * (( *nb_eigen_all) + 1)];
      }
      F(j, j) = (*alpha)*(fact.values())[j + j * (( *nb_eigen_all) + 1)]
			     + beta[j] * (fact.values())[(j + 1) + j * (( *nb_eigen_all) + 1)];

      if(j + 1 <= *nb_eigen){
        F(j + 1, j) = beta[j] * (fact.values())[(j + 1) + (j + 1) * (( *nb_eigen_all) + 1)];
      }
    }

  /* set the vectors*/

  Teuchos::SerialDenseVector<int, double> rhs((*nb_eigen) + 1);

  /*set the solution to zero*/
  rhs = 0.0;
  /* rhs[0] must be setted to beta*/
  rhs(0) = (fact.values())[0];

  /*create the lsqr solver context and set it up*/

  int info2;
  int rank;
  char TRANS = 'N';
  int lwork;
  double* work;
  double wkopt;

  lapack.GELS(TRANS,3,3,1,F.values(),F.stride(),rhs.values(), rhs.stride(), &wkopt, -1, &info2);
  lwork = (int)wkopt;
  work = new double [lwork];
  lapack.GELS(TRANS,3,3,1,F.values(),F.stride(),rhs.values(), rhs.stride(), work, lwork, &info2);

  //rhs sequence vector store the solution of least squares problem
  for(i = 0; i < *nb_eigen_all; i++){
    eta[i] = rhs.values[i];
  }

  //free the allocate memory
  for(i = 0; i < (*nb_eigen_all) + 3; i++){
    delete [] gamma[i];
  }

  delete [] gamma;

  for(i = 0; i < (*nb_eigen_all) + 1; i++){
    delete [] mm_tmp[i];
  }

  delete [] mm_tmp;

  delete [] work;

  return 0;
}
