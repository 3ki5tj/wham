#ifndef MBAR_H__
#define MBAR_H__



/* multistate Bennett's acceptance ratio method
 * this module requires `xvg.h` */



#include <time.h>
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
  const xdouble *beta; /* temperature array, reference */
  xdouble *tot; /* total number of visits to each temperature */
  xdouble *lntot; /* logarithm of the total number of visits to each temperature */
  xdouble *res; /* difference between new and old lnz */
  xdouble **lnden; /* logarithm of the denominator */
  xvg_t **xvg; /* xvg trajectory */
  unsigned flags;
} mbar_t;



/* flags */
#define MBAR_RMCOM    0x0010 /* normalize by removing the center of mass motion */
#define MBAR_NOEST    0x0020 /* do not estimate lnz at the beginning */



static mbar_t *mbar_open(int nbeta, const xdouble *beta,
    xvg_t **xvg, unsigned flags)
{
  mbar_t *mbar;
  int i;

  xnew(mbar, 1);
  mbar->nbeta = nbeta;
  mbar->beta = beta;
  mbar->xvg = xvg;
  xnew(mbar->tot, nbeta);
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



static void mbar_close(mbar_t *mbar)
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
__inline static xdouble mbar_lnadd(xdouble a, xdouble b)
{
  xdouble c;

  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + LOG(1 + EXP(-c));
}



/* shift `lnz` such that `lnz[0] == 0` */
static void mbar_shift(xdouble *lnz, int nbeta)
{
  int i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



/* normalize by remove the center of mass motion */
static void mbar_rmcom(xdouble *arr, const xdouble *tot, int n)
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
static void mbar_normalize(xdouble *lnz, int nbeta, const void *ptr)
{
  const mbar_t *w = (const mbar_t *) ptr;

  if ( w->flags & MBAR_RMCOM ) {
    mbar_rmcom(lnz, w->tot, nbeta);
  } else {
    mbar_shift(lnz, nbeta);
  }
}



/* estimate the partition function using the single histogram method */
static void mbar_estimatelnz(mbar_t *mbar, xdouble *lnz)
{
  xvg_t *xvg;
  int i, j, n, nbeta = mbar->nbeta;
  xdouble db, e, dlnz;

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
    lnz[j] = lnz[j - 1] + (n > 0 ? LOG(n) - dlnz : 0);
  }

  mbar_normalize(lnz, nbeta, mbar);
}



static xdouble mbar_step(mbar_t *mbar, xdouble *lnz, xdouble *res, xdouble damp)
{
  xvg_t *xvg;
  int i, j, l, n, nbeta = mbar->nbeta;
  xdouble lny, lnden, e, err = 0;

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

  if ( damp > 0 ) {
    for ( i = 0; i < nbeta; i++ ) {
      lnz[i] += damp * res[i];
    }

    mbar_normalize(lnz, nbeta, mbar);
  }

  return err;
}



/* weighted histogram analysis method */
static xdouble mbar(int nbeta,
    xvg_t **xvg, const xdouble *beta, xdouble *lnz, unsigned flags,
    xdouble damp, int itmin, int itmax, xdouble tol, int verbose)
{
  mbar_t *mbar = mbar_open(nbeta, beta, xvg, flags);
  int it;
  xdouble err, errp;
  clock_t t0, t1;

  if ( !(flags & MBAR_NOEST) ) {
    mbar_estimatelnz(mbar, lnz);
  }

  t0 = clock();
  err = errp = 1e30;
  for ( it = 0; it < itmax; it++ ) {
    err = mbar_step(mbar, lnz, mbar->res, damp);
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
  fprintf(stderr, "MBAR converged in %d steps, error %g, time %.4fs\n",
      it, (double) err, 1.0*(t1 - t0)/CLOCKS_PER_SEC);
  mbar_close(mbar);
  return err;
}



#ifdef ENABLE_MDIIS
/* MDIIS method */
#include "../mdiis_xdbl.h"



static xdouble mbar_getres(void *mbar, xdouble *lnz, xdouble *res)
{
  return mbar_step((mbar_t *) mbar, lnz, res, 0);
}



static xdouble mbar_mdiis(int nbeta,
    xvg_t **xvg, const xdouble *beta, xdouble *lnz, unsigned flags,
    int nbases, xdouble damp, int queue, xdouble threshold,
    int itmin, int itmax, xdouble tol, int verbose)
{
  mbar_t *mbar = mbar_open(nbeta, beta, xvg, flags);
  xdouble err;

  if ( !(flags & MBAR_NOEST) ) {
    mbar_estimatelnz(mbar, lnz);
  }

  err = iter_mdiis(lnz, nbeta,
      mbar_getres, mbar_normalize, mbar,
      nbases, damp, queue, threshold,
      itmin, itmax, tol, NULL, verbose);
  mbar_close(mbar);
  return err;
}

#endif /* ENABLE_MDIIS */



static xdouble mbarx(int nbeta, xvg_t **xvg,
    const xdouble *beta, xdouble *lnz, unsigned flags,
    xdouble damp, int nbases, int update_method, xdouble threshold,
    int itmin, int itmax, xdouble tol, int verbose, int method)
{
  if ( method == MBAR_DIRECT ) {
    return mbar(nbeta, xvg, beta, lnz, flags,
        damp, itmin, itmax, tol, verbose);
#ifdef ENABLE_MDIIS
  } else if ( method == MBAR_MDIIS ) {
    return mbar_mdiis(nbeta, xvg, beta, lnz, flags,
        nbases, damp, update_method, threshold,
        itmin, itmax, tol, verbose);
#endif
  }

  return 0;
}



#endif /* MBAR_H__ */

