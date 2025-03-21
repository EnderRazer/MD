#ifndef RANDOM_NORMALIZED_H
#define RANDOM_NORMALIZED_H

#include <cmath>

inline double ran2(int *idum) {
  const int IM1 = 2147483563;
  const int IM2 = 2147483399;
  const double AM = (1.0 / IM1);
  const int IMM1 = (IM1 - 1);
  const int IA1 = 40014;
  const int IA2 = 40692;
  const int IQ1 = 53668;
  const int IQ2 = 52774;
  const int IR1 = 12211;
  const int IR2 = 3791;
  const int NTAB = 32;
  const double NDIV = (1 + double(IMM1) / double(NTAB));
  const double EPS1 = 1.2e-7;
  const double RNMX = (1.0 - EPS1);

  int j;
  int k;
  static int idum2 = 123456789;
  static int iy = 0;
  static int iv[NTAB];
  double tempran;

  if (*idum <= 0) { /* initialize */
    if (-(*idum) < 1) {
      *idum = 1; /* be sure to prevent idum = 0 */
    } else {
      *idum = -(*idum);
    }
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) { /*load shuffle table (after 8 warm-ups) */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k * IQ1) - k * IR1;
      if (*idum < 0) {
        *idum += IM1;
      }
      if (j < NTAB) {
        iv[j] = *idum;
      }
    }
    iy = iv[0];
  }
  k = (*idum) / IQ1; /* start here where not init*/
  *idum = IA1 * (*idum - k * IQ1) - k * IR1;
  if (*idum < 0) {
    *idum += IM1;
  }
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0) {
    idum2 += IM2;
  }
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) {
    iy += IMM1;
  }
  if ((tempran = AM * iy) > RNMX) {
    return RNMX;
  } else
    return tempran;
}

inline double gasdev(int idum) {
  int *idum2 = &idum;
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if (iset == 0) {
    do {
      v1 = 2.0 * ran2(idum2) - 1.0;
      v2 = 2.0 * ran2(idum2) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}

#endif
