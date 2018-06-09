#ifndef _CONVHULL_H_
#define _CONVHULL_H_

#include <complex>
#include <algorithm>

int convhull(std::complex<double> *ab, std::complex<double> *c, std::complex<double> *d, int n, int *ne, int offset, int mu){

  std::complex<double> *hk = new std::complex<double>[n + 2];
  double s, cs, dzero, deps;
  int m, i, j, m1;
  int *l = new int[n + 2];

  /* Offset:  Acces aux elements [1+offset:n+offset] des tableaux a et b */
  /* mu: 	! Acces aux elements [1+mu:n+mu] des tableaux cr, ci, dr, di */

  dzero = 0.0;
  deps = 1.0e-16;

  for(i = 0; i < n + 2; i++){
    l[i] = 1;
  }

  /* trouve l'indice de la valeur propres à plus petite partie réelle*/
  s = ab[offset].real();
  m1 = 0;

  for(i = 1; i < n; i++){
    if(ab[i + offset].real() <= s){
      s = ab[offset].real();
      m1 = i;
    }
  }

  *ne = 0;
  hk[0].real(s); hk[0].imag(dzero); /* init complex array by s+i*0 */

  /* met la partie imaginaire de la valeur propre de plus petite partie réelle à 0 */
  if(ab[m1 + offset].imag() < deps){
    l[m1] = 0;
  } else{
    ab[m1 + offset].real(s);
    ab[m1 + offset].imag(dzero);
  }

  for(i = 1; i <= n; ++i){
    s = -1.0;
    m = m1;
    for(j = 0; j < n; j++){
      if(l[j]){
        cs = sqrt(std::pow(ab[j + offset].real() - hk[i - 1].real(),2) +
                  std::pow(ab[j + offset].imag() - hk[i - 1].imag(),2));
        if(cs >= deps){
          cs = ab[j + offset].imag() - ( hk[i - 1].imag() ) / cs;
        } else {
          continue;
        }
        if(cs - s > -deps){
          if(cs - s < deps){
            if(ab[j + offset].real() > ab[m + offset].real() + deps){
              m = j;
            }
            if(s >= (1 - deps) && ab[j + offset].imag() > ab[m + offset].imag()){
              m = j;
            }
          } else {
            s = cs;
            m = j;
          }
        }
      }
      if( m == m1){
        break;
      } else{
        hk[i] = ab[m + offset];
        l[m] = 0;
        *ne = i;
        c[i - 1 + mu].real((hk[i].real() + hk[i - 1].real()) / 2);
        c[i - 1 + mu].imag((hk[i].imag() + hk[i - 1].imag()) / 2);

        d[i - 1 + mu].real((hk[i].real() - hk[i - 1].real()) / 2);
        d[i - 1 + mu].imag((hk[i].imag() - hk[i - 1].imag()) / 2);

        for(j = 0; j < n; j++){
          if(l[j]){
            if(ab[j + offset].real() <= hk[i].real()){
              l[j] = 0;
            }
          }
        }
       }
    }
  }

  if(hk[*ne].imag() != dzero){
    c[(*ne) + mu].real(hk[*ne].real());
    c[(*ne) + mu].imag(hk[*ne].imag() / 2);

    d[(*ne) + mu].real(dzero);
    d[(*ne) + mu].imag(hk[*ne].imag() / 2);

    *ne = (*ne) + 1;
  }

  delete [] l;
  delete [] hk;

  return 0;
}


#endif
