#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Version.hpp"


int main(int argc, char* argv[])
{

  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

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
