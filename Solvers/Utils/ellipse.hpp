#ifndef _ELLIPSE_H_
#define _ELLIPSE_H_

#include <complex>
#include <algorithm>
#include <cmath>


int ellipse3(std::complex<double> xy1, std::complex<double> xy2, std::complex<double> xy3, double * a2, double * b2, double * c, int * info){

  double d2, z, r;

  r = ( real(xy2) - real(xy1) ) * (real(xy3) - real(xy2)) * (real(xy3) - real(xy1));

  if(r > 1.0e-14){
    z = std::pow(imag(xy1), 2) * (real(xy2) - real(xy3)) +
        std::pow(imag(xy2), 2) * (real(xy3) - real(xy1)) +
        std::pow(imag(xy3), 2) * (real(xy1) - real(xy2));

        if(std::fabs(z) < 1.0e-14){
          return 1;
        }

        *c = ( (std::pow(imag(xy1),2)) * ( (std::pow(real(xy2),2)) - (std::pow(real(xy3),2))) +
               (std::pow(imag(xy2),2)) * ( (std::pow(real(xy3),2)) - (std::pow(real(xy1),2))) +
               (std::pow(imag(xy3),2)) * ( (std::pow(real(xy1),2)) - (std::pow(real(xy2),2))) ) / (2 * z);

        *a2 = std::pow(*c, 2) - (std::pow(imag(xy1),2)*real(xy2)*real(xy3)*(real(xy2)-real(xy3))
			     +std::pow(imag(xy2),2)*real(xy1)*real(xy3)*(real(xy3)-real(xy1))
			     +std::pow(imag(xy3),2)*real(xy1)*real(xy2)*(real(xy1)-real(xy2)))/z;

        d2 = (*a2) * (1.0 - z / r);

        *b2 = *a2 - d2;
  } else {
    return 1;
  }

  *info = 1;
  if(*a2 < 0.0 || *b2 < 0.0){
    *info = 0;
  }

  return 0;
}

int ellipse2(std::complex<double> xy1, std::complex<double> xy2, double * a2, double * b2, double * c, int * info){
	double a, b, s, t, z, q, d2, sign;

	a = (real(xy2) - real(xy1))/2;
	b = (real(xy2) + real(xy1))/2;
	s = (imag(xy2) - imag(xy1))/2;
	t = (imag(xy2) + imag(xy1))/2;

	if(t <= 1.e-15 || std::fabs(a) <= 1.e-15)	{
	  return 1;
	}

	if(std::fabs(s) >= 1.e-15){
		q = (s / t+t / s) / 2;
		sign = 1.0;

		if(a*s <= 0) {
      sign = -1.0;
    }

		z = ((-q) + sign * std::sqrt(std::pow(q,2) + 3)) * (a / 3);
		if(std::fabs(z)<1.e-15){
		  return 1;
		}
		*a2 = (z + a * t / s) * (z + a * s / t);
		d2 = (*a2) * (z - s * t / a) / z;
		*b2 = (*a2) - d2;
		*c = b + z;
	} else {
		*a2 = 2 * std::pow(a,2);
		*b2 = 2 * std::pow(*c,2);
		d2 = (*a2) - (*b2);
		z = 0.0;
		*c = b + z;
	}

	*info = 1;
	if((*a2) < 0.0 || (*b2) < 0.0) {
    	*info = 0;
  	}
	return 0;
}

int test(std::complex<double> * hk, int n, double a2, double b2, double c, int * info){
	int i;
	double d;

	*info = 1;
	for(i = 0; i < n; i++){
		d = b2 * std::pow(real(hk[i])-c,2) + a2 * std::pow(imag(hk[i]),2) - a2 * b2;
		if(d >= 1.e-6){
			*info = 0;
			return 1;
		}
	}
	return 0;
}

void ellipse(std::complex<double> * c, std::complex<double> * d, int n, int mu, double * co, double * ao2, double * do2, double * dr, int * info){
	int i, j, k, info1, info2;
	std::complex<double> * hk  = new std::complex<double> [n];
	double bo2, a2, b2, cc, aire;
	printf("info===%d\n", *info);

	std::cout << "$} Ellipse Allocating work memory " << n << "\n";

    std::cout << "$} Ellipse Computing Edges\n";
	/* CALCUL DES SOMMETS A PARTIR DES CENTRES ET DES DEMI-DISTANCES */
	for(i = 0;i < n - 1; i++){
		hk[i] = c[i] - d[i];
		printf("hk[%d] = %f+%fi\n",i,real(hk[i]), imag(hk[i]));
	}

	printf("----------------------\n");
	if(mu != 0 && mu < n){
		hk[n-1] = c[mu] + d[mu];
		printf("hk[%d] = %f+%fi\n",n-1,real(hk[n-1]), imag(hk[n-1]));
	}

	hk[n-1] = c[n-2] + d[n-2];

    std::cout << "$} Ellipse Research two points optimal ellipse\n";

	/* Recherche de l'ellipse optimale a deux points */
	info1 = 0;
	info2 = 0;

	*ao2 = 0.0;
	bo2 = 0.0;
	*co = 0.0;

	aire = 1.e16;

	for(i = 0;i < n - 1; i++){
		for(j = i + 1;j < n; j++){
			info1 = 0;
			ellipse2(hk[i],hk[j],&a2,&b2,&cc,&info1);
            if(info1 == 0) continue;
			test(hk,n,a2,b2,cc,&info1);
			if(info1 == 1){
				if(a2 * b2 < aire){
					info2 = 1;
					(*ao2) = a2;
					bo2 = b2;
					(*co) = cc;
					aire = (*ao2) * bo2;
					printf("aire = %f\n", aire);
				}
			}
		}
	}

    std::cout << "$} Ellipse Research three points optimal ellipse\n";

	/* Recherche de l'ellipse optimale a trois points */
	if(info2 != 1){
		aire = 1.e16;

		for(i = 0; i < n - 2; i++){
			for(j = i + 1; j < n - 1; j++){
				for(k = j + 1; k < n; k++){
					*info = 0;
					ellipse3(hk[i],hk[j],hk[k],&a2,&b2,&cc,info);
					if(*info == 0) continue;
					test(hk,n,a2,b2,cc,info);
					if(*info == 1){
						if(a2 * b2 < aire){
							*ao2 = a2;
							bo2 = b2;
							*co = cc;
							aire = (*ao2) * b2;
						}
					}
				}
			}
		}
	}

	*do2 = (*ao2) - bo2;
	*dr= (*ao2) >= bo2;
	*ao2 = std::sqrt(*ao2);

	printf("*info = = %d\n", *info);

	delete [] hk;

}


#endif
