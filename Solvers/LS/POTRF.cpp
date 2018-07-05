#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Version.hpp"
#include "../Utils/libs.hpp"
#include "../Utils/convhull.hpp"
#include "../Utils/ellipse.hpp"
#include <iostream>
#include "precond.hpp"
#include "lsp.hpp"


int main(int argc, char* argv[])
{

    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

    std::cout << "Hello, World!" << std::endl;

    double eps = epsilon();

    std::cout << "eps = " << eps << std::endl;

    int length = 5;
    int length2 = 5;
    std::complex<double> *arr = new std::complex<double>[length];
    std::complex<double> *arr2 = new std::complex<double>[length2];

    int i;

    arr2[0].real(-0.830999);
    arr2[0].imag(0.514128);

    arr2[1].real(-0.774165);
    arr2[1].imag(0.424414);

    arr2[2].real(-0.463353);
    arr2[2].imag(0.363388);

    arr2[3].real(-0.444916);
    arr2[3].imag(0.517988);

    arr2[4].real(-0.287419);
    arr2[4].imag(0.356902);


    for( i = 0; i < length2; i++){
      std::cout << "arr2[" << i << "] = " << arr2[i] << '\n';
    }

    std::cout << "===== tri operation =====" << '\n';

    int ch_signe;
    tri(arr2, length2, &ch_signe);
    for( i = 0; i < length2; i++){
      std::cout << "arr2[" << i << "] = " << arr2[i] << '\n';
    }

    std::cout << "ch_signe = " << ch_signe << '\n';
   
    std::complex<double> * c = new std::complex<double>[20];
    std::complex<double> * d = new std::complex<double>[20];

    int mu1;

    convhull(arr2, c, d, ch_signe, &mu1, 0, 0);

    std::cout << "mu1 = " << mu1 << '\n';

/*
    for(i = 0; i < 20; i++){
        printf("c[%d] = %f + %fi, d[%d] = %f + %fi\n", i, c[i].real(), c[i].imag(),i, d[i].real(), d[i].imag());
    }
*/

    double a_ell, c_ell, d_ell, d_reel;
    int info3;

    ellipse(c, d, mu1 + 1, mu1, &c_ell, &a_ell, &d_ell, &d_reel, &info3);
    printf("info3 = %d, c_ell = %f, a_ell = %f, d_ell = %f, d_reel = %f \n", info3, c_ell, a_ell, d_ell,d_reel);

    //LSPrecond(a_ell, d_ell,c_ell,eta, &alpha, beta, delta, c, d,&mu, &ls_eigen, &ls_eigen_min, &eigen_max);


  // Creating an instance of the LAPACK class for double-precision routines looks like:
  Teuchos::LAPACK<int, double> lapack;

  // This instance provides the access to all the LAPACK routines.

  Teuchos::SerialDenseMatrix<int, double> My_Matrix(3,3);
  Teuchos::SerialDenseVector<int, double> My_Vector(3);

  My_Vector.random();

  // Print out the original linear system.
  std::cout << "ORIGINAL MATRIX:" << std::endl;

  //row 0
  My_Matrix(0,0) = 4.0; My_Matrix(0,1) = 12.0; My_Matrix(0,2) = -16.0;

  My_Matrix(1,0) = 12.0; My_Matrix(1,1) = 37.0; My_Matrix(1,2) = -43.0;

  My_Matrix(2,0) = -16.0; My_Matrix(2,1) = -43.0; My_Matrix(2,2) = 98;

  std::cout << My_Matrix << std::endl;

  int info;
  char UPLO = 'L';
  lapack.POTRF(UPLO,3,My_Matrix.values(),My_Matrix.stride(),&info );

  // Print out the solution.
  std::cout << "CHOLESKY FACTORIZED MATRIX FROM LAPACK:" << std::endl;

  std::cout << My_Matrix << std::endl;

/*
  Teuchos in Trilinos 12.12.1

  ORIGINAL MATRIX:


  Values_copied : yes
  Rows : 3
  Columns : 3
  LDA : 3
  4 12 -16
  12 37 -43
  -16 -43 98

  CHOLESKY FACTORIZED MATRIX FROM LAPACK:


  Values_copied : yes
  Rows : 3
  Columns : 3
  LDA : 3
  2 12 -16
  6 1 -43
  -8 5 3

  The upper part of my_matrix is its cholesky facotrization
  */

  /*Solve a least square problem*/
  
  std::cout << "The given right hand side is :" << std::endl;

  std::cout << My_Vector << std::endl;

  int info2;
  int rank;
  char TRANS = 'N';
  int lwork;
  double* work;
  double wkopt;

  lapack.GELS(TRANS,3,3,1,My_Matrix.values(),My_Matrix.stride(),My_Vector.values(), My_Vector.stride(), &wkopt, -1, &info2);
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  lapack.GELS(TRANS,3,3,1,My_Matrix.values(),My_Matrix.stride(),My_Vector.values(), My_Vector.stride(), work, lwork, &info2);

  std::cout << My_Vector << std::endl;


  return 0;
}
