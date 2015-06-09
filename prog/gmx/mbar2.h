#ifndef MBAR_H__
#define MBAR_H__



/* multistate Bennett's acceptance ratio method
 * this module requires `xvg.h` */



#include <time.h>
#include <float.h>
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
  const double *bx, *by; /* temperature array, reference */
  double *res; /* difference between new and old lnz */
  double *lntot; /* total number of visits to each temperature */
  double **lnden; /* logarithm of the denominator */
  xvg_t **xvg; /* xvg trajectory */
} mbar2_t;



static mbar2_t *mbar2_open(int nbeta,
    const double *bx, const double *by, xvg_t **xvg)
{
  mbar2_t *mbar;
  int i;

  xnew(mbar, 1);
  mbar->nbeta = nbeta;
  mbar->bx = bx;
  mbar->by = by;
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



static void mbar2_close(mbar2_t *mbar)
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
__inline static double mbar2_lnadd(double a, double b)
{
  double c;

  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + log(1 + exp(-c));
}



/* estimate the partition function using the single histogram method */
static void mbar2_estimatelnz(mbar2_t *mbar, double *lnz)
{
  xvg_t *xvg;
  int i, k, kk, l, n, nbeta = mbar->nbeta;
  double db, dbx, dby, x, y, dlnz;

  lnz[0] = 0;
  for ( k = 1; k < nbeta; k++ ) {
    db = DBL_MAX;
    kk = 0;
    /* select the nearest neighbor of k as the reference */
    for ( l = 0; l < k; l++ ) {
      if ( l == k ) continue;
      dbx = mbar->bx[k] - mbar->bx[l];
      dby = mbar->by[k] - mbar->by[l];
      if ( (x = dbx * dbx + dby * dby) < db ) {
        db = x;
        kk = l;
      }
    }
    dbx = mbar->bx[k] - mbar->bx[kk];
    dby = mbar->by[k] - mbar->by[kk];

    dlnz = LOG0;
    /* Z(kk) / Z(k)
     * = < exp( [bx(k) - bx(kk)] x + [by(k) - by(kk)] y ) >_k
     **/
    xvg = mbar->xvg[k];
    n = xvg->n;
    for ( i = 0; i < n; i++ ) {
      x = xvg->y[0][i];
      y = xvg->y[2][i];
      dlnz = mbar2_lnadd(dlnz, dbx * x + dby * y);
    }
    lnz[k] = lnz[kk] + (n > 0 ? log(n) - dlnz : 0);
  }
}



static void mbar2_normalize(double *lnz, int nbeta)
{
  int i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



static double mbar2_step(mbar2_t *mbar, double *lnz, double *res, int update)
{
  xvg_t *xvg;
  int i, j, l, n, nbeta = mbar->nbeta;
  double lny, lnden, x, y, err;

  /* compute the denominator of each frame */
  for ( i = 0; i < nbeta; i++ ) {
    xvg = mbar->xvg[i];
    n = xvg->n;
    for ( l = 0; l < n; l++ ) {
      /* compute the denominator for frame (i, l) */
      x = xvg->y[0][l];
      y = xvg->y[2][l];
      lnden = LOG0;
      for ( j = 0; j < nbeta; j++ ) {
        lnden = mbar2_lnadd(lnden,
            mbar->lntot[j] - mbar->bx[j] * x - mbar->by[j] * y - lnz[j]);
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
        x = xvg->y[0][l];
        y = xvg->y[2][l];
        lny = mbar2_lnadd(lny,
            -mbar->bx[i] * x - mbar->by[j] * y - mbar->lnden[j][l]);
      }
    }
    res[i] = lny - lnz[i];
  }

  mbar2_normalize(res, nbeta);
  for ( err = 0, i = 0; i < nbeta; i++ ) {
    if ( fabs(res[i]) > err ) {
      err = fabs(res[i]);
    }
  }

  if ( update ) {
    for ( i = 0; i < nbeta; i++ ) {
      lnz[i] += res[i];
    }
    mbar2_normalize(lnz, nbeta);
  }

  return err;
}



/* weighted histogram analysis method */
static double mbar2(int nbeta,
    xvg_t **xvg, const double *bx, const double *by, double *lnz,
    int itmax, double tol, int verbose)
{
  mbar2_t *mbar = mbar2_open(nbeta, bx, by, xvg);
  int it;
  double err, errp;
  clock_t t0, t1;

  mbar2_estimatelnz(mbar, lnz);

  t0 = clock();
  err = errp = 1e30;
  for ( it = 0; it < itmax; it++ ) {
    err = mbar2_step(mbar, lnz, mbar->res, 1);
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
  fprintf(stderr, "MBAR2 converged in %d steps, error %g, time %.4fs\n",
      it, err, 1.0*(t1 - t0)/CLOCKS_PER_SEC);
  mbar2_close(mbar);
  return err;
}



#ifdef ENABLE_MDIIS
/* MDIIS method */
#include "../mdiis.h"



static double mbar2_getres(void *mbar, double *lnz, double *res)
{
  return mbar2_step((mbar2_t *) mbar, lnz, res, 0);
}



static double mbar2_mdiis(int nbp, xvg_t **xvg,
    const double *bx, const double *by, double *lnz,
    int nbases, double damp, int queue, double threshold,
    int itmax, double tol, int verbose)
{
  mbar2_t *mbar = mbar2_open(nbp, bx, by, xvg);
  double err;

  mbar2_estimatelnz(mbar, lnz);
  err = iter_mdiis(lnz, nbp,
      mbar2_getres, mbar2_normalize, mbar,
      nbases, damp, queue, threshold, itmax, tol, verbose);
  mbar2_close(mbar);
  return err;
}

#endif /* ENABLE_MDIIS */



#endif /* MBAR_H__ */

