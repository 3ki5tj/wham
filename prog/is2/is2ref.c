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



static double is2_exactx(int lx, int ly, double beta,
    double *eav, double *var, double *skw)
{
  double lnz, db = 1e-6;
  double beta1, eav1, var1;
  double beta2, eav2, var2;

  lnz = is2_exact(lx, ly, beta, eav, var);
  *var /= beta * beta;
  beta1 = beta - db;
  is2_exact(lx, ly, beta1, &eav1, &var1);
  var1 /= beta1 * beta1;
  beta2 = beta + db;
  is2_exact(lx, ly, beta2, &eav2, &var2);
  var2 /= beta2 * beta2;
  *skw = (var1 - var2) / (2 * db);
  return lnz;
}



/* Gaussian partition method (not very successful)
 * seek the position where the energy distributions
 * at two temperatures are the same */
static double getdlnzgp(double beta1, double eav1, double var1,
    double beta2, double eav2, double var2)
{
  double a1 = 1/var1, a2 = 1/var2;
  double A, B, C, D, num;

  A = a2 - a1;
  B = a2 * eav2 - a1 * eav1;
  /* the term log(a1/a2) is better dropped
   * for thermodyanmic consistency */
  C = a2 * eav2 * eav2 - a1 * eav1 * eav1; /* + log(a1/a2); */
  D = sqrt(B*B - A*C);
  if ( eav1 > eav2 ) num = B + D;
  else num = B - D;
  return num/A * (beta1 - beta2);
}



/* Tilted Gaussian density of states */
static double getdlnztg(double beta1, double eav1, double var1,
    double beta2, double eav2, double var2)
{
  double de = eav2 - eav1, db = beta2 - beta1;
  return -(eav1 + eav2) / 2 * db + (1/var2 - 1/var1) * de * de / 12;
}



/* ln var approximation */
static double getdlnzlnv(double beta1, double eav1, double var1,
    double beta2, double eav2, double var2)
{
  double de = eav2 - eav1, db = beta2 - beta1;
  return -(eav1 + eav2) / 2 * db + log(var2 / var1) * de * db / 12;
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
  double *beta, *eav, *var, *skw, *lnz;
  double *lnza, *lnzb, *lnzc, *lnzuip, *lnzui, *lnzgp, *lnztg, *lnzlnv;
  double db, dlnza, dlnzb, dlnzc, dlnzuip, dlnzgp, dlnztg, dlnzlnv;
  int iT, jT;

  model_default_is2(m);
  model_doargs(m, argc, argv);

  xnew(beta, m->nT);
  xnew(eav, m->nT);
  xnew(var, m->nT);
  xnew(skw, m->nT);
  xnew(lnz, m->nT);
  xnew(lnza, m->nT);
  xnew(lnzb, m->nT);
  xnew(lnzc, m->nT);
  xnew(lnzuip, m->nT);
  xnew(lnzui, m->nT);
  xnew(lnzgp, m->nT);
  xnew(lnztg, m->nT);
  xnew(lnzlnv, m->nT);

  /* set up the temperatures */
  for ( iT = 0; iT < m->nT; iT++ ) {
    double T = m->Tmin + m->Tdel * iT;
    beta[iT] = 1 / T;
    lnz[iT] = is2_exactx(m->L, m->L, beta[iT], &eav[iT], &var[iT], &skw[iT]);
  }
  for ( iT = m->nT - 1; iT >= 0; iT-- ) {
    lnz[iT] -= lnz[0];
  }

  /* compute pairwise estimates values */
  lnza[0] = lnzb[0] = lnzc[0] = lnzuip[0] = 0;
  lnzgp[0] = lnztg[0] = lnzlnv[0] = 0;
  for ( iT = 1; iT < m->nT; iT++ ) {
    jT = iT - 1;
    db = beta[iT] - beta[jT];
    dlnza = -(eav[iT] + eav[jT]) / 2 * db;
    lnza[iT] = lnza[jT] + dlnza;
    dlnzb = dlnza - (var[iT] - var[jT]) / 12 * db * db;
    lnzb[iT] = lnzb[jT] + dlnzb;
    dlnzc = dlnza + (skw[iT] + skw[jT]) / 24 * db * db * db;
    lnzc[iT] = lnzc[jT] + dlnzc;
    dlnzuip = getlnzui(2, lnzui, beta + jT, eav + jT, var + jT);
    lnzuip[iT] = lnzuip[jT] + dlnzuip;
    dlnzgp = getdlnzgp(beta[jT], eav[jT], var[jT], beta[iT], eav[iT], var[iT]);
    lnzgp[iT] = lnzgp[jT] + dlnzgp;
    dlnztg = getdlnztg(beta[jT], eav[jT], var[jT], beta[iT], eav[iT], var[iT]);
    lnztg[iT] = lnztg[jT] + dlnztg;
    dlnzlnv = getdlnzlnv(beta[jT], eav[jT], var[jT], beta[iT], eav[iT], var[iT]);
    lnzlnv[iT] = lnzlnv[jT] + dlnzlnv;
  }

  getlnzui(m->nT, lnzui, beta, eav, var);
  for ( iT = 0; iT < m->nT; iT++ ) {
    printf("%5.3f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f\n",
        1/beta[iT], lnz[iT], lnza[iT], lnzb[iT], lnzc[iT],
        lnzuip[iT], lnzui[iT], lnzgp[iT], lnztg[iT], lnzlnv[iT]);
  }

  free(beta);
  free(eav);
  free(var);
  free(skw);
  free(lnz);
  free(lnza);
  free(lnzb);
  free(lnzc);
  free(lnzuip);
  free(lnzui);
  free(lnzgp);
  free(lnztg);
  free(lnzlnv);

  return 0;
}

