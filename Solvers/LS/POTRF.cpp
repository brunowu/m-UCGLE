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

    int eigen_max = 10, ls_eigen = 10, ls_eigen_min = 10;

    double *delta = new double [eigen_max + 1];
    double *eta = new double [eigen_max + 1];
    double *beta = new double [eigen_max + 1];

    double alpha;

    int mu;


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

    double a_ell, c_ell, d_ell, d_reel;
    int info3;

    ellipse(c, d, mu1 + 1, mu1, &c_ell, &a_ell, &d_ell, &d_reel, &info3);
    printf("info3 = %d, c_ell = %f, a_ell = %f, d_ell = %f, d_reel = %f \n", info3, c_ell, a_ell, d_ell,d_reel);

    mu = mu1;


    std::cout << "mu = " << mu << '\n';

    LSPrecond(a_ell, d_ell,c_ell,eta, &alpha, beta, delta, c, d,&mu, &ls_eigen, &ls_eigen_min, &eigen_max);

    for(i = 0; i < eigen_max; i++){
      printf("eta[%d] = %f, alpha = %f, beta[%d] = %f, deta[%d] =  %f\n", i, eta[i], alpha, i, beta[i],i, delta[i]);
    }

  return 0;
}
