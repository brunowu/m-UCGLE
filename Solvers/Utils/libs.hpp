#ifndef _LIBS_H_
#define _LIBS_H_

#include <complex>
#include <algorithm>
#include <cmath>


/* calcul de l'epsilon machine */
double epsilon(void){
	double one = 1.0, two = 2.0,temp = 1.0;

	while(one+temp>one){
    temp/=two;
  }
	return two*temp;
}

/*compare func for c++ sort*/
bool compare(std::complex<double> a, std::complex<double> b){

	if(a.real() == b.real()){
		return a.imag() < b.imag();
	}

	return real(a) < real(b);

}

/*sort of array*/
void sort(std::complex<double> *arr, int beg, int end){

	std::sort(arr + beg, arr + end, compare);

}

/*sort the eigenvalues in ascending order depending on the real parts
Find the rank of the last one with a negative real part*/
void tri(std::complex<double> *vp, int nValues, int *ch_signe){
	int i;
	/*sort the eigenvalues depending on the realparts*/
	for(i = 0; i < nValues; ++i){
		sort(vp, 0, nValues);

		*ch_signe = 0;

		while(vp[*ch_signe].real() < 0.0){
			(*ch_signe)++;
			if(*ch_signe > nValues){
				break;
			}
		}
	}
}

/*remove the eigenvalues whose imaginairy part is negetif*/
int epurer(std::complex<double> *vp, int *nValues){
	int i, nKeep = 0;
	double tmp;
	printf("EPURE FUNC with *nVlaues = %d\n",*nValues);

	if(*nValues <= 0){return 1;}
	std::complex<double> *eigen_keep = new std::complex<double> [*nValues];

	for(i = 0; i < *nValues; i++){
		if(vp[i].real() < 0.0){
			tmp = -vp[i].real();
			//printf("tmp = %f, i = %d\n", tmp, i);
		}
		else{
			tmp = vp[i].real();
		}

		if(vp[i].imag() >= 0.0 && tmp != 0.0 && tmp > 1e-06){
			eigen_keep[nKeep].real(vp[i].real());
			eigen_keep[nKeep].imag(vp[i].imag());

			nKeep++;
		}
	}

	for(i = 0; i < nKeep; i++){
		vp[i].real(eigen_keep[i].real());
		vp[i].imag(eigen_keep[i].imag());
	}

	*nValues = nKeep;

	delete [] eigen_keep;

	return 0;
}
#endif
