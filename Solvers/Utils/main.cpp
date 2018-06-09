#include <iostream>
#include "libs.hpp"
#include "convhull.hpp"
#include "ellipse.hpp"

int main()
{
    std::cout << "Hello, World!" << std::endl;

    double eps = epsilon();

    std::cout << "eps = " << eps << std::endl;

    int length = 5;
    std::complex<double> *arr = new std::complex<double>[length];
    int i;
    for( i = 0; i < length; i++){
      arr[i].real((i - 1)^2);
      arr[i].imag((i - 1)^2);
      std::cout << "arr[" << i << "] = " << arr[i] << '\n';
    }

    std::cout << "===== tri operation =====" << '\n';

    int ch_signe;
    tri(arr, length, &ch_signe);
    for( i = 0; i < length; i++){
      std::cout << "arr[" << i << "] = " << arr[i] << '\n';
    }

    std::cout << "ch_signe = " << ch_signe << '\n';

    std::cout << "===== epurer operation =====" << '\n';

    epurer(arr, &length);
    std::cout << "new length = " << length << '\n';
    for( i = 0; i < length; i++){
      std::cout << "new arr[" << i << "] = " << arr[i] << '\n';
    }

    return 0;
}
