#ifndef WHAM2_H__
#define WHAM2_H__



/* Two-dimensional weighted histogram analysis method
 * this module requires `hist2.h` */



#include <time.h>
#include "hist2.h"



#ifndef LOG0
#define LOG0 -1e9
#endif



enum { WHAM_DIRECT = 0, WHAM_MDIIS = 1, WHAM_NMETHODS };
const char *wham_methods[] = {"Direct", "MDIIS"};



typedef struct {
  const double *bx, *by; /* temperature array, reference */
  double *tot; /* total number of visits to each temperature */
  double *lntot; /* logarithm of the total number of visits to each temperature */
  double *res; /* difference between new and old lnz */
  double *htot; /* overall histogram */
  double *lndos; /* density of states */
  hist2_t *hist; /* histograms, reference */
  int imin, imax, jmin, jmax;
  unsigned flags;
} wham2_t;



/* flags */
#define WHAM2_RMCOM   0x0010 /* normalize by removing the center of mass motion */
#define WHAM2_NOEST   0x0020 /* do not estimate lnz at the beginning */



static wham2_t *wham2_open(const double *bx, const double *by,
    hist2_t *hist, unsigned flags)
{
  wham2_t *w;
  int i, j, ij, n = hist->n, m = hist->m, nm = n * m, k, nbeta = hist->rows;

  xnew(w, 1);
  w->bx = bx;
  w->by = by;
  w->hist = hist;
  xnew(w->tot, hist->rows);
  xnew(w->lntot, hist->rows);
  xnew(w->res, hist->rows);
  xnew(w->htot, nm);
  xnew(w->lndos, nm);

  /* compute the total */
  for ( ij = 0; ij < nm; ij++ ) {
    w->htot[ij] = 0;
  }
  for ( k = 0; k < nbeta; k++ ) {
    double x, s = 0;
    for ( ij = 0; ij < nm; ij++ ) {
      x = hist->arr[k * nm + ij];
      s += x;
      w->htot[ij] += x;
    }
    w->tot[k] = s;
    w->lntot[k] = (s > 0) ? log(s) : LOG0;
  }

  /* determine the boundaries */
  /* find imin */
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < m; j++ ) {
      if ( w->htot[i * m + j] > 0 ) break;
    }
    if ( j < m ) break;
  }
  w->imin = i;

  /* find imax */
  for ( i = n - 1; i >= 0; i-- ) {
    for ( j = 0; j < m; j++ ) {
      if ( w->htot[i * m + j] > 0 ) break;
    }
    if ( j < m ) break;
  }
  w->imax = i + 1;

  /* find jmin */
  for ( j = 0; j < m; j++ ) {
    for ( i = 0; i < n; i++ ) {
      if ( w->htot[i * m + j] > 0 ) break;
    }
    if ( i < n ) break;
  }
  w->jmin = j;

  /* find jmax */
  for ( j = m - 1; j >= 0; j-- ) {
    for ( i = 0; i < n; i++ ) {
      if ( w->htot[i * m + j] > 0 ) break;
    }
    if ( i < n ) break;
  }
  w->jmax = j + 1;

  w->flags = flags;

  fprintf(stderr, "i: [%d, %d), j: [%d, %d)\n", w->imin, w->imax, w->jmin, w->jmax);
  return w;
}



static void wham2_close(wham2_t *w)
{
  free(w->tot);
  free(w->lntot);
  free(w->res);
  free(w->htot);
  free(w->lndos);
  free(w);
}



/* log(exp(a) + exp(b)) */
__inline static double wham2_lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + log(1 + exp(-c));
}



/* save the density of states to file `fn` */
static int wham2_savelndos(wham2_t *w, const char *fn)
{
  hist2_t *hist = w->hist;
  FILE *fp;
  int i, j, ij, m = hist->m;
  double x, y;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  for ( i = w->imin; i < w->imax; i++ ) {
    x = hist->xmin + (i + .5) * hist->dx;
    for ( j = w->jmin; j < w->jmax; j++ ) {
      y = hist->ymin + (j + .5) * hist->dy;
      ij = i*m + j;
      if ( w->lndos[ij] > LOG0 )
        fprintf(fp, "%g %g %g\n", x, y, w->lndos[ij]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}



/* compute the average properties from the density of states */
static void wham2_getav(wham2_t *w, const char *fn)
{
  hist2_t *hist = w->hist;
  double bx, by, x, y, xx, xy, yy, lnw, lnx, lny;
  double lnz, slnx, slny, slnxx, slnxy, slnyy;
  int i, j, ij, m = hist->m, nbeta = hist->rows, k;
  FILE *fp = NULL;

  if ( fn != NULL && (fp = fopen(fn, "w")) == NULL )
    fprintf(stderr, "cannot write %s\n", fn);

  for ( k = 0; k < nbeta; k++ ) {
    bx = w->bx[k];
    by = w->by[k];
    lnz = slnx = slny = slnxx = slnxy = slnyy = LOG0;
    for ( i = w->imin; i < w->imax; i++ ) {
      x = (i + .5) * hist->dx;
      lnx = log(x);
      for ( j = w->jmin; j < w->jmax; j++ ) {
        ij = i*m + j;
        if ( w->lndos[ij] <= LOG0 ) continue;
        y = (j + .5) * hist->dy;
        lny = log(y);

        lnw = w->lndos[ij] - bx * x - by * y;
        lnz = wham2_lnadd(lnz, lnw);
        slnx = wham2_lnadd(slnx, lnx + lnw);
        slny = wham2_lnadd(slny, lny + lnw);
        slnxx = wham2_lnadd(slnxx, lnx + lnx + lnw);
        slnxy = wham2_lnadd(slnxy, lnx + lny + lnw);
        slnyy = wham2_lnadd(slnyy, lny + lny + lnw);
      }
    }
    x = exp(slnx - lnz);
    y = exp(slny - lnz);
    xx = exp(slnxx - lnz) - x * x;
    xy = exp(slnxy - lnz) - x * y;
    yy = exp(slnyy - lnz) - y * y;
    x += hist->xmin;
    y += hist->ymin;
    if ( fp != NULL )
      fprintf(fp, "%g %g %g %g %g %g %g %g\n",
          bx, by, x, y, xx, xy, yy,
          lnz - bx * hist->xmin - by * hist->ymin);
  }
  if ( fp != NULL ) fclose(fp);
}



/* shift `lnz` such that `lnz[0] == 0` */
static void wham2_shift(double *lnz, int nbeta)
{
  int i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



/* normalize by removing the center of mass motion */
static void wham2_rmcom(double *arr,
    const double *tot, int n)
{
  int i;
  double s = 0, sy = 0;

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
static void wham2_normalize(double *lnz, int nbeta, const void *ptr)
{
  const wham2_t *w = (const wham2_t *) ptr;

  if ( w->flags & WHAM2_RMCOM ) {
    wham2_rmcom(lnz, w->tot, nbeta);
  } else {
    wham2_shift(lnz, nbeta);
  }
}



/* estimate the partition function using the single histogram method */
static void wham2_estimatelnz(wham2_t *w, double *lnz)
{
  hist2_t *hist = w->hist;
  int i, j, ij, n = hist->n, m = hist->m, nm = n * m;
  int nbeta = hist->rows, k, l, kk;
  double dbx, dby, db, x, y, h, s, dlnz;

  lnz[0] = 0;
  for ( k = 1; k < nbeta; k++ ) {
    db = DBL_MAX;
    kk = 0;
    /* select the nearest neighbor of k as the reference */
    for ( l = 0; l < k; l++ ) {
      if ( l == k ) continue;
      dbx = w->bx[k] - w->bx[l];
      dby = w->by[k] - w->by[l];
      if ( (x = dbx * dbx + dby * dby) < db ) {
        db = x;
        kk = l;
      }
    }
    dbx = w->bx[k] - w->bx[kk];
    dby = w->by[k] - w->by[kk];

    s = 0;
    dlnz = LOG0;
    /* Z(kk) / Z(k)
     * = < exp( [bx(k) - bx(kk)] x + [by(k) - by(kk)] y ) >_k
     *    Sum_{x, y} h_k(x, y) exp( [bx(k) - by(kk)] x + [by(k) - by(kk)] y )
     * = ---------------------------------------------------------------------
     *                     Sum_{x, y} h_k(x, y)
     **/
    for ( i = w->imin; i < w->imax; i++ ) {
      x = hist->xmin + (i + .5) * hist->dx;
      for ( j = w->jmin; j < w->jmax; j++ ) {
        ij = i*m + j;

        h = hist->arr[k*nm + ij];
        if ( h <= 0 ) continue;
        y = hist->ymin + (j + .5) * hist->dy;

        s += h;
        dlnz = wham2_lnadd(dlnz, log(h) + dbx * x + dby * y);
      }
    }
    lnz[k] = lnz[kk] + (s > 0 ? log(s) - dlnz : 0);
  }

  wham2_normalize(lnz, nbeta, w);
}



/* compute the partition function from the density of states */
static void wham2_getlnz(wham2_t *w, double *lnz)
{
  hist2_t *hist = w->hist;
  int i, j, ij, k, m = hist->m, nbeta = hist->rows;
  double bx, by, x, y;

  for ( k = 0; k < nbeta; k++ ) {
    bx = w->bx[k];
    by = w->by[k];
    lnz[k] = LOG0;
    for ( i = w->imin; i < w->imax; i++ ) {
      x = hist->xmin + (i + .5) * hist->dx;
      for ( j = w->jmin; j < w->jmax; j++ ) {
        ij = i*m + j;
        if ( w->lndos[ij] <= LOG0) continue;
        y = hist->ymin + (j + .5) * hist->dy;

        /* res[k] is the new partition function */
        lnz[k] = wham2_lnadd(lnz[k], w->lndos[ij] - bx * x - by * y);
      }
    }
  }

  wham2_normalize(lnz, nbeta, w);
}



static double wham2_step(wham2_t *w, double *lnz, double *res,
    double damp)
{
  hist2_t *hist = w->hist;
  int i, j, ij, ijmin, n = hist->n, m = hist->m, nm = n * m;
  int k, nbeta = hist->rows;
  double x, y, bx, by, lnden, err;

  ijmin = -1;
  for ( i = w->imin; i < w->imax; i++ ) {
    x = hist->xmin + (i + .5) * hist->dx;
    for ( j = w->jmin; j < w->jmax; j++ ) {
      y = hist->ymin + (j + .5) * hist->dy;
      ij = i*m + j;
      if ( w->htot[ij] <= 0 ) continue;

      lnden = LOG0;
      /*        num           Sum_k h_k(x, y)
       * dos = ----- = ---------------------------------------------
       *        den     Sum_k tot_k exp(-bx_k * x - by_k * y) / Z_k
       * */
      for ( k = 0; k < nbeta; k++ ) {
        bx = w->bx[k];
        by = w->by[k];
        lnden = wham2_lnadd(lnden, w->lntot[k] - bx * x - by * y - lnz[k]);
      }
      w->lndos[ij] = log(w->htot[ij]) - lnden;
    }
  }

  /* shift the baseline of the density of states */
  for ( x = w->lndos[ijmin], ij = 0; ij < nm; ij++ )
    if ( w->lndos[ij] > LOG0 )
      w->lndos[ij] -= x;

  /* refresh the partition function */
  wham2_getlnz(w, res);

  for ( err = 0, i = 0; i < nbeta; i++ ) {
    res[i] -= lnz[i];
    if ( fabs(res[i]) > err ) {
      err = fabs(res[i]);
    }
  }

  if ( damp > 0 ) {
    for ( i = 0; i < nbeta; i++ ) {
      lnz[i] += damp * res[i];
    }
    wham2_normalize(lnz, nbeta, w);
  }

  return err;
}



/* iteratively compute the logarithm of the density of states
 * using the weighted histogram method */
static double wham2_getlndos(wham2_t *w, double *lnz,
    double damp, int itmin, int itmax, double tol, int verbose)
{
  int it;
  double err, errp;
  clock_t t0, t1;

  t0 = clock();
  err = errp = 1e30;
  for ( it = 1; it <= itmax; it++ ) {
    err = wham2_step(w, lnz, w->res, damp);
    if ( verbose ) {
      fprintf(stderr, "it %d, err %g -> %g\n",
          it, errp, err);
    }
    if ( err < tol && it > itmin ) {
      break;
    }
    errp = err;
  }

  t1 = clock();
  fprintf(stderr, "WHAM2 converged in %d steps, error %g, time %.4fs\n",
      it, err, 1.0*(t1 - t0)/CLOCKS_PER_SEC);
  return err;
}



/* two-dimensional weighted histogram analysis method */
static double wham2(hist2_t *hist,
    const double *bx, const double *by, double *lnz,
    unsigned flags, double damp,
    int itmin, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham2_t *w = wham2_open(bx, by, hist, flags);
  double err;

  if ( !(flags & WHAM2_NOEST) ) {
    wham2_estimatelnz(w, lnz);
  }

  err = wham2_getlndos(w, lnz,
      damp, itmin, itmax, tol, verbose);
  if ( fnlndos ) {
    wham2_savelndos(w, fnlndos);
  }
  wham2_getav(w, fneav);
  wham2_close(w);
  return err;
}



#ifdef ENABLE_MDIIS
/* MDIIS method */
#include "mdiis.h"



static double wham2_getres(void *w, double *lnz, double *res)
{
  return wham2_step((wham2_t *) w, lnz, res, 0);
}



static double wham2_mdiis(hist2_t *hist,
    const double *bx, const double *by, double *lnz,
    unsigned flags,
    int nbases, double damp, int queue, double threshold,
    int itmin, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham2_t *w = wham2_open(bx, by, hist, flags);
  double err;

  if ( !(flags & WHAM2_NOEST) ) {
    wham2_estimatelnz(w, lnz);
  }

  err = iter_mdiis(lnz, hist->rows,
      wham2_getres, wham2_normalize, w,
      nbases, damp, queue, threshold,
      itmin, itmax, tol, NULL, verbose);
  if ( fnlndos ) wham2_savelndos(w, fnlndos);
  wham2_getav(w, fneav);
  wham2_close(w);
  return err;
}

#endif /* ENABLE_MDIIS */



static double wham2x(hist2_t *hist,
    const double *bx, const double *by, double *lnz,
    unsigned flags,
    double damp, int nbases, int update_method, double threshold,
    int itmin, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav, int method)
{
  if ( method == WHAM_DIRECT ) {
    return wham2(hist, bx, by, lnz, flags,
        damp, itmin, itmax, tol,
        verbose, fnlndos, fneav);
#ifdef ENABLE_MDIIS
  } else if ( method == WHAM_MDIIS ) {
    return wham2_mdiis(hist, bx, by, lnz, flags,
        nbases, damp, update_method, threshold,
        itmin, itmax, tol, verbose,
        fnlndos, fneav);
#endif
  }

  return 0;
}



#endif /* WHAM2_H__ */

