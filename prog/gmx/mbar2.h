#ifndef MBAR2_H__
#define MBAR2_H__



/* multistate Bennett's acceptance ratio method
 * this module requires `xvg.h` */



#include <time.h>
#include <float.h>
#include "xvg.h"
#include "../xdouble.h"



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



enum { MBAR_DIRECT = 0, MBAR_MDIIS = 1, MBAR_NMETHODS };
const char *mbar_methods[] = {"Direct", "MDIIS"};



typedef struct {
  int nbeta;
  const xdouble *bx, *by; /* temperature array, reference */
  xdouble *tot; /* total number of visits to each temperature */
  xdouble *lntot; /* logarithm of the total number of visits to each temperature */
  xdouble *res; /* difference between new and old lnz */
  xdouble **lnden; /* logarithm of the denominator */
  xvg_t **xvg; /* xvg trajectory */
  unsigned flags;
} mbar2_t;



/* flags */
#define MBAR2_RMCOM   0x0010 /* normalize by removing the center of mass motion */
#define MBAR2_NOEST   0x0020 /* do not estimate lnz at the beginning */



static mbar2_t *mbar2_open(int nbeta,
    const xdouble *bx, const xdouble *by,
    xvg_t **xvg, unsigned flags)
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
    xdouble x = xvg[i]->n;
    mbar->tot[i] = x;
    mbar->lntot[i] = (x > 0) ? LOG(x) : LOG0;
  }

  mbar->flags = flags;

  return mbar;
}



static void mbar2_close(mbar2_t *mbar)
{
  int i;

  for ( i = 0; i < mbar->nbeta; i++ ) {
    free(mbar->lnden[i]);
  }
  free(mbar->tot);
  free(mbar->lntot);
  free(mbar->lnden);
  free(mbar->res);
  free(mbar);
}



/* log(exp(a) + exp(b)) */
__inline static xdouble mbar2_lnadd(xdouble a, xdouble b)
{
  xdouble c;

  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + LOG(1 + EXP(-c));
}



/* shift `lnz` such that `lnz[0] == 0` */
static void mbar2_shift(xdouble *lnz, int nbeta)
{
  int i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



/* normalize by remove the center of mass motion */
static void mbar2_rmcom(xdouble *arr, const xdouble *tot, int n)
{
  int i;
  xdouble s = 0, sy = 0;

  for ( i = 0; i < n; i++ ) {
    s += tot[i];
    sy += arr[i] * tot[i];
  }
  sy /= s;
  for ( i = 0; i < n; i++ ) {
    arr[i] -= sy;
  }
}



/* normalize lnz */
static void mbar2_normalize(xdouble *lnz, int nbeta, const void *ptr)
{
  const mbar2_t *w = (const mbar2_t *) ptr;

  if ( w->flags & MBAR2_RMCOM ) {
    mbar2_rmcom(lnz, w->tot, nbeta);
  } else {
    mbar2_shift(lnz, nbeta);
  }
}



/* estimate the partition function using the single histogram method */
static void mbar2_estimatelnz(mbar2_t *mbar, xdouble *lnz)
{
  xvg_t *xvg;
  int i, k, kk, l, n, nbeta = mbar->nbeta;
  xdouble db, dbx, dby, x, y, dlnz;

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
      y = xvg->y[1][i];
      dlnz = mbar2_lnadd(dlnz, dbx * x + dby * y);
    }
    lnz[k] = lnz[kk] + (n > 0 ? LOG(n) - dlnz : 0);
  }

  mbar2_normalize(lnz, nbeta, mbar);
}



static xdouble mbar2_step(mbar2_t *mbar, xdouble *lnz, xdouble *res, xdouble damp)
{
  xvg_t *xvg;
  int i, j, l, n, nbeta = mbar->nbeta;
  xdouble lny, lnden, x, y, err;

  /* compute the denominator of each frame */
  for ( i = 0; i < nbeta; i++ ) {
    xvg = mbar->xvg[i];
    n = xvg->n;
    for ( l = 0; l < n; l++ ) {
      /* compute the denominator for frame (i, l) */
      x = xvg->y[0][l];
      y = xvg->y[1][l];
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
        y = xvg->y[1][l];
        lny = mbar2_lnadd(lny,
            -mbar->bx[i] * x - mbar->by[j] * y - mbar->lnden[j][l]);
      }
    }
    res[i] = lny - lnz[i];
  }

  mbar2_normalize(res, nbeta, mbar);
  for ( err = 0, i = 0; i < nbeta; i++ ) {
    if ( fabs(res[i]) > err ) {
      err = fabs(res[i]);
    }
  }

  if ( damp > 0 ) {
    for ( i = 0; i < nbeta; i++ ) {
      lnz[i] += damp * res[i];
    }
    mbar2_normalize(lnz, nbeta, mbar);
  }

  return err;
}



/* weighted histogram analysis method */
static xdouble mbar2(int nbeta,
    xvg_t **xvg, const xdouble *bx, const xdouble *by, xdouble *lnz,
    unsigned flags,
    xdouble damp, int itmin, int itmax, xdouble tol, int verbose)
{
  mbar2_t *mbar = mbar2_open(nbeta, bx, by, xvg, flags);
  int it;
  xdouble err, errp;
  clock_t t0, t1;

  if ( !(flags & MBAR2_NOEST) ) {
    mbar2_estimatelnz(mbar, lnz);
  }

  t0 = clock();
  err = errp = 1e30;
  for ( it = 0; it < itmax; it++ ) {
    err = mbar2_step(mbar, lnz, mbar->res, damp);
    if ( verbose ) {
      fprintf(stderr, "it %d, err %g -> %g\n",
          it, (double) errp, (double) err);
    }
    if ( err < tol && it > itmin ) {
      break;
    }
    errp = err;
  }
  t1 = clock();
  fprintf(stderr, "MBAR2 converged in %d steps, error %g, time %.4fs\n",
      it, (double) err, 1.0*(t1 - t0)/CLOCKS_PER_SEC);
  mbar2_close(mbar);
  return err;
}



#ifdef ENABLE_MDIIS
/* MDIIS method */
#include "../mdiis_xdbl.h"



static xdouble mbar2_getres(void *mbar, xdouble *lnz, xdouble *res)
{
  return mbar2_step((mbar2_t *) mbar, lnz, res, 0);
}



static xdouble mbar2_mdiis(int nbp, xvg_t **xvg,
    const xdouble *bx, const xdouble *by, xdouble *lnz, unsigned flags,
    int nbases, xdouble damp, int update_method, xdouble threshold,
    int itmin, int itmax, xdouble tol, int verbose)
{
  mbar2_t *mbar = mbar2_open(nbp, bx, by, xvg, flags);
  xdouble err;

  if ( !(flags & MBAR2_NOEST) ) {
    mbar2_estimatelnz(mbar, lnz);
  }

  err = iter_mdiis(lnz, nbp,
      mbar2_getres, mbar2_normalize, mbar,
      nbases, damp, update_method, threshold,
      itmin, itmax, tol, NULL, verbose);
  mbar2_close(mbar);
  return err;
}

#endif /* ENABLE_MDIIS */



static xdouble mbar2x(int nbp, xvg_t **xvg,
    const xdouble *bx, const xdouble *by, xdouble *lnz, unsigned flags,
    xdouble damp, int nbases, int update_method, xdouble threshold,
    int itmin, int itmax, xdouble tol, int verbose, int method)
{
  if ( method == MBAR_DIRECT ) {
    return mbar2(nbp, xvg, bx, by, lnz, flags,
        damp, itmin, itmax, tol, verbose);
#ifdef ENABLE_MDIIS
  } else if ( method == MBAR_MDIIS ) {
    return mbar2_mdiis(nbp, xvg, bx, by, lnz, flags,
        nbases, damp, update_method, threshold,
        itmin, itmax, tol, verbose);
#endif
  }

  return 0;
}



#endif /* MBAR2_H__ */

