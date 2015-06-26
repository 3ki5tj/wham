/* reference values for two-dimensional Ising model */
#define IS2_LB 6
#include "is2.h"
#include <time.h>
#define IS2_MODEL
#include "../whammodel.h"



#ifndef xnew
#define xnew(x, n) { int i_; \
  if ((x = calloc((size_t) (n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %ld\n", #x, (long) n); \
    exit(1); \
  } \
  for ( i_ = 0; i_ < n; i_ ++ ) x[i_] = 0; }
#endif



static void model_default_is2(model_t *m)
{
  model_default(m);
  m->nT = 1600;
  m->Tmin = 1.5;
  m->Tdel = 0.001;
  m->L = 1 << IS2_LB;
}



/* compute lnz from umbrella integration */
static double getlnzui(int m, double *lnz,
    double *beta, double *eav, double *var)
{
  double de, emin, emax;
  double *logn, *bst, ymax, num, den, du, x, dx;
  double *barr, *sarr;
  int i, j, ne;

  /* estimate the range and width of energy histograms */
  emin = emax = eav[0];
  de = dx = var[0];
  for ( j = 1; j < m; j++ ) {
    if ( var[j] < de ) de = var[j];
    if ( var[j] > dx ) dx = var[j];
    if ( eav[j] < emin ) emin = eav[j];
    if ( eav[j] > emax ) emax = eav[j];
    //printf("j %d, beta %g, eav %g\n", j, beta[j], eav[j]);
  }
  de = sqrt(de);
  dx = sqrt(dx);
  emin -= 10 * de;
  emax += 10 * dx;
  de = de / 500;
  ne = (int) ((emax - emin) / de);
  //printf("de %g, emin %g, emax %g\n", de, emin, emax);

  xnew(barr, ne);
  xnew(sarr, ne);
  xnew(logn, m);
  xnew(bst, m);

  for ( i = 0; i < ne; i++ ) {
    x = emin + i * de;
    ymax = -1e300;
    for ( j = 0; j < m; j++ ) {
      dx = x - eav[j];
      logn[j] = -dx * dx / 2 / var[j] - 0.5 * log(2 * M_PI * var[j]);
      if ( logn[j] > ymax ) {
        ymax = logn[j];
      }
      bst[j] = beta[j] - dx / var[j];
    }
    num = den = 0;
    for ( j = 0; j < m; j++ ) {
      du = exp(logn[j] - ymax);
      num += bst[j] * du;
      den += du;
    }
    barr[i] = num/den;
  }

  /* integrate the statistical temperature */
  sarr[0] = 0;
  for ( i = 1; i < ne; i++ ) {
    sarr[i] = sarr[i-1] + (barr[i] + barr[i-1]) * de/2;
  }

  /* compute the partition function */
  for ( j = 0; j < m; j++ ) {
    lnz[j] = -1e30;
  }
  for ( i = 0; i < ne; i++ ) {
    x = emin + i * de;
    for ( j = 0; j < m; j++ ) {
      lnz[j] = lnadd(lnz[j], sarr[i] - beta[j] * x);
    }
  }
  /* shift the origin of the partition function */
  for ( j = m - 1; j >= 0; j-- ) {
    lnz[j] -= lnz[0];
  }

  free(barr);
  free(sarr);
  free(logn);
  free(bst);

  return lnz[m-1] - lnz[0];
}



int main(int argc, char **argv)
{
  model_t m[1];
  double *beta, *eav, *var, *lnz;
  double *lnza, *lnzb, *lnzuip, *lnzui;
  double db, dlnza, dlnzb, dlnzuip;
  int iT, jT;

  model_default_is2(m);
  model_doargs(m, argc, argv);

  xnew(beta, m->nT);
  xnew(eav, m->nT);
  xnew(var, m->nT);
  xnew(lnz, m->nT);
  xnew(lnza, m->nT);
  xnew(lnzb, m->nT);
  xnew(lnzuip, m->nT);
  xnew(lnzui, m->nT);

  /* set up the temperatures */
  for ( iT = 0; iT < m->nT; iT++ ) {
    double T = m->Tmin + m->Tdel * iT;
    beta[iT] = 1 / T;
    lnz[iT] = is2_exact(m->L, m->L, beta[iT], &eav[iT], &var[iT]);
    var[iT] *= T * T;
  }
  for ( iT = m->nT - 1; iT >= 0; iT-- ) {
    lnz[iT] -= lnz[0];
  }

  /* compute pairwise estimates values */
  lnza[0] = lnzb[0] = lnzuip[0];
  for ( iT = 1; iT < m->nT; iT++ ) {
    jT = iT - 1;
    db = beta[iT] - beta[jT];
    dlnza = -(eav[iT] + eav[jT]) / 2 * db;
    lnza[iT] = lnza[jT] + dlnza;
    dlnzb = dlnza - (var[iT] - var[jT]) / 12 * db * db;
    lnzb[iT] = lnzb[jT] + dlnzb;
    dlnzuip = getlnzui(2, lnzui, beta + jT, eav + jT, var + jT);
    lnzuip[iT] = lnzuip[jT] + dlnzuip;
  }

  getlnzui(m->nT, lnzui, beta, eav, var);
  for ( iT = 0; iT < m->nT; iT++ ) {
    printf("%5.3f %14.7f %14.7f %14.7f %14.7f %14.7f\n",
        1/beta[iT], lnz[iT], lnza[iT], lnzb[iT], lnzuip[iT], lnzui[iT]);
  }

  free(beta);
  free(eav);
  free(var);
  free(lnz);
  free(lnza);
  free(lnzb);
  free(lnzuip);
  free(lnzui);

  return 0;
}

