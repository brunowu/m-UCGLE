#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#define PI 3.14159265

using namespace std;

int main () {

        cout << "]> Starting to generate the eigenvalues\n";

        ofstream myfile ("spectrum.txt");
        int N = 200,i;

        srand((unsigned)time(0));

        double random_value, random_value2;

        if (myfile.is_open()){
                myfile << "%%MatrixMarket matrix coordinate real general\n";
                myfile << N << " " << N << " "<< N << "\n";
                for (i=0;i<N;i++){
                        random_value =  (double)rand()/RAND_MAX*0.001+0.00001;
                        random_value2 = (double)rand()/RAND_MAX*0.001;
                        myfile << i+1 << " "<< random_value << " "<<  random_value2 << "\n";
                }
        cout << "]> Finsh to generate the eigenvalues\n";
        }
        else cout << "Unable to open file";

        return 0;
}

