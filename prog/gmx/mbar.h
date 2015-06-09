#ifndef MBAR_H__
#define MBAR_H__



/* multistate Bennett's acceptance ratio method
 * this module requires `xvg.h` */



#include <time.h>
#include "xvg.h"



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((size_t) (n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %ld\n", #x, (long) n); \
    exit(1); \
  }}
#endif



#ifndef LOG0
#define LOG0 -1e9
#endif



typedef struct {
  int nbeta;
  const double *beta; /* temperature array, reference */
  double *res; /* difference between new and old lnz */
  double *lntot; /* total number of visits to each temperature */
  double **lnden; /* logarithm of the denominator */
  xvg_t **xvg; /* xvg trajectory */
} mbar_t;



static mbar_t *mbar_open(int nbeta, const double *beta, xvg_t **xvg)
{
  mbar_t *mbar;
  int i;

  xnew(mbar, 1);
  mbar->nbeta = nbeta;
  mbar->beta = beta;
  mbar->xvg = xvg;
  xnew(mbar->lntot, nbeta);
  xnew(mbar->res, nbeta);
  xnew(mbar->lnden, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    xnew(mbar->lnden[i], xvg[i]->n);
  }

  /* copy the total */
  for ( i = 0; i < nbeta; i++ ) {
    double x = xvg[i]->n;
    mbar->lntot[i] = (x > 0) ? log(x) : LOG0;
  }

  return mbar;
}



static void mbar_close(mbar_t *mbar)
{
  int i;

  for ( i = 0; i < mbar->nbeta; i++ ) {
    free(mbar->lnden[i]);
  }
  free(mbar->lnden);
  free(mbar->res);
  free(mbar->lntot);
  free(mbar);
}



/* log(exp(a) + exp(b)) */
__inline static double mbar_lnadd(double a, double b)
{
  double c;

  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + log(1 + exp(-c));
}



/* estimate the partition function using the single histogram method */
static void mbar_estimatelnz(mbar_t *mbar, double *lnz)
{
  xvg_t *xvg;
  int i, j, n, nbeta = mbar->nbeta;
  double db, e, dlnz;

  lnz[0] = 0;
  for ( j = 1; j < nbeta; j++ ) {
    /* estimate the free energy different between
     * the pair j - 1 and j */
    db = mbar->beta[j] - mbar->beta[j - 1];
    dlnz = LOG0;
    /* Z(j-1) / Z(j) = < exp( [beta(j) - beta(j-1)] E ) >_j */
    xvg = mbar->xvg[j];
    n = xvg->n;
    for ( i = 0; i < n; i++ ) {
      e = xvg->y[0][i];
      dlnz = mbar_lnadd(dlnz, db * e);
    }
    lnz[j] = lnz[j - 1] + (n > 0 ? log(n) - dlnz : 0);
  }
}



static void mbar_normalize(double *lnz, int nbeta)
{
  int i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



static double mbar_step(mbar_t *mbar, double *lnz, double *res, int update)
{
  xvg_t *xvg;
  int i, j, l, n, nbeta = mbar->nbeta;
  double lny, lnden, e, err = 0;

  /* compute the denominator of each frame */
  for ( i = 0; i < nbeta; i++ ) {
    xvg = mbar->xvg[i];
    n = xvg->n;
    for ( l = 0; l < n; l++ ) {
      /* compute the denominator for frame (i, l) */
      e = xvg->y[0][l];
      lnden = LOG0;
      for ( j = 0; j < nbeta; j++ ) {
        lnden = mbar_lnadd(lnden, mbar->lntot[j] - mbar->beta[j] * e - lnz[j]);
      }
      mbar->lnden[i][l] = lnden;
    }
  }

  for ( i = 0; i < nbeta; i++ ) {
    /* compute the residue of lnz[i] */
    lny = LOG0;
    for ( j = 0; j < nbeta; j++ ) {
      xvg = mbar->xvg[j];
      n = xvg->n;
      for ( l = 0; l < n; l++ ) {
        e = xvg->y[0][l];
        lny = mbar_lnadd(lny, -mbar->beta[i] * e - mbar->lnden[j][l]);
      }
    }
    res[i] = lny - lnz[i];
    if ( fabs(res[i]) > err ) {
      err = fabs(res[i]);
    }
  }

  if ( update ) {
    for ( i = 0; i < nbeta; i++ ) {
      lnz[i] += res[i];
    }
    mbar_normalize(lnz, nbeta);
  }

  return err;
}



/* weighted histogram analysis method */
static double mbar(int nbeta,
    xvg_t **xvg, const double *beta, double *lnz,
    int itmax, double tol, int verbose)
{
  mbar_t *mbar = mbar_open(nbeta, beta, xvg);
  int it;
  double err, errp;
  clock_t t0, t1;

  mbar_estimatelnz(mbar, lnz);

  t0 = clock();
  err = errp = 1e30;
  for ( it = 0; it < itmax; it++ ) {
    err = mbar_step(mbar, lnz, mbar->res, 1);
    if ( verbose ) {
      fprintf(stderr, "it %d, err %g -> %g\n",
          it, errp, err);
    }
    if ( err < tol ) {
      break;
    }
    errp = err;
  }
  t1 = clock();
  fprintf(stderr, "MBAR converged in %d steps, error %g, time %.4fs\n",
      it, err, 1.0*(t1 - t0)/CLOCKS_PER_SEC);
  mbar_close(mbar);
  return err;
}



#ifdef ENABLE_MDIIS
/* MDIIS method */
#include "../mdiis.h"



static double mbar_getres(void *mbar, double *lnz, double *res)
{
  return mbar_step((mbar_t *) mbar, lnz, res, 0);
}



static double mbar_mdiis(int nbeta,
    xvg_t **xvg, const double *beta, double *lnz,
    int nbases, double damp, int queue, double threshold,
    int itmax, double tol, int verbose)
{
  mbar_t *mbar = mbar_open(nbeta, beta, xvg);
  double err;

  mbar_estimatelnz(mbar, lnz);
  err = iter_mdiis(lnz, nbeta,
      mbar_getres, mbar_normalize, mbar,
      nbases, damp, queue, threshold, itmax, tol, verbose);
  mbar_close(mbar);
  return err;
}

#endif /* ENABLE_MDIIS */



#endif /* MBAR_H__ */

