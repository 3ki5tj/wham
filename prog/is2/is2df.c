/* comparison of approximate free energy formulas */
#include "is2.h"
#include <stdlib.h>
#include <math.h>



static double getdfui(double beta1, double eav1, double var1,
    double beta2, double eav2, double var2)
{
  double de, de1, de2, emin, emax;
  double logn1, logn2, du, x, dx1, dx2, b1, b2;
  double *barr, *sarr, lnz1, lnz2;
  int i, ne;

  de1 = sqrt(var1);
  de2 = sqrt(var2);
  de = (de1 + de2)/2;
  emin = eav2 - 10 * de;
  emax = eav1 + 10 * de;
  de = de / 200;
  ne = (int) ((emax - emin) / de);
  barr = calloc(ne, sizeof(*barr));
  sarr = calloc(ne, sizeof(*sarr));

  for ( i = 0; i < ne; i++ ) {
    x = emin + i * de;
    dx1 = x - eav1;
    dx2 = x - eav2;
    logn1 = -dx1 * dx1 / var1 - 0.5 * log(2*M_PI*var1);
    logn2 = -dx2 * dx2 / var2 - 0.5 * log(2*M_PI*var2);
    b1 = -dx1 / var1 + beta1;
    b2 = -dx2 / var2 + beta2;
    if ( logn1 > logn2 ) {
      du = exp(logn2 - logn1);
      barr[i] = (b1 + b2 * du) / (1 + du);
    } else {
      du = exp(logn1 - logn2);
      barr[i] = (b1 * du + b2) / (du + 1);
    }
  }

  /* integrate the statistical temperature */
  sarr[0] = 0;
  for ( i = 1; i < ne; i++ ) {
    sarr[i] = sarr[i-1] + (barr[i] + barr[i-1]) * de/2;
  }

  /* compute the free energy difference */
  lnz1 = lnz2 = -1e30;
  for ( i = 0; i < ne; i++ ) {
    x = emin + i * de;
    lnz1 = lnadd(lnz1, sarr[i] - beta1 * x);
    lnz2 = lnadd(lnz2, sarr[i] - beta2 * x);
  }

  free(barr);
  free(sarr);

  return lnz1 - lnz2;
}



int main(int argc, char **argv)
{
  double T1 = 2.0, beta1, lnz1, eav1, var1;
  double beta2, lnz2, eav2, var2;
  double db, c, x, y, ich2, df, df0, dfser, dfui;
  int l = 64, i;

  if ( argc >= 2 ) {
    T1 = atof(argv[1]);
  }
  if ( argc >= 3 ) {
    l = atoi(argv[2]);
  }

  beta1 = 1 / T1;
  lnz1 = is2_exact(l, l, beta1, &eav1, &var1);
  var1 /= beta1 * beta1;

  for ( i = 1; i <= 20; i++ ) {
    db = 0.01 * i;
    beta2 = beta1 + db;
    lnz2 = is2_exact(l, l, beta2, &eav2, &var2);
    var2 /= beta2 * beta2;
    c = (eav1 - eav2) / 2;
    x = db * c / 2;
    y = exp(-x);
    ich2 = 2 * y / (1 + y * y);
    ich2 *= ich2;
    df = -(lnz2 - lnz1);
    df0 = (eav2 + eav1) / 2 * db;
    dfser = df0 + (var2 - var1) * db * db / 12;
    dfui = getdfui(beta1, eav1, var1, beta2, eav2, var2);
    printf("%5.3f %12.7f %12.7f %12.7f %12.7f %g\n",
        db, df, df0, dfser, dfui, x*2);
  }

  return 0;
}

