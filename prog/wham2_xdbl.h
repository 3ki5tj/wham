#ifndef WHAM2_H__
#define WHAM2_H__



/* Two-dimensional weighted histogram analysis method
 * this module requires `hist2.h` */



#include <time.h>
#include "hist2.h"
#include "xdouble.h"



#ifndef LOG0
#define LOG0 -1e9
#endif



enum { WHAM_DIRECT = 0, WHAM_MDIIS = 1, WHAM_NMETHODS };
const char *wham_methods[] = {"Direct", "MDIIS"};



typedef struct {
  const xdouble *bx, *by; /* temperature array, reference */
  xdouble *tot; /* total number of visits to each temperature */
  xdouble *lntot; /* logarithm of the total number of visits to each temperature */
  xdouble *res; /* difference between new and old lnz */
  xdouble *htot; /* overall histogram */
  xdouble *lndos; /* density of states */
  hist2_t *hist; /* histograms, reference */
  int imin, imax, jmin, jmax;
  unsigned flags;
} wham2_t;



/* flags */
#define WHAM2_RMCOM   0x0010 /* normalize by removing the center of mass motion */
#define WHAM2_NOEST   0x0020 /* do not estimate lnz at the beginning */



static wham2_t *wham2_open(const xdouble *bx, const xdouble *by,
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
    xdouble x, s = 0;
    for ( ij = 0; ij < nm; ij++ ) {
      x = hist->arr[k * nm + ij];
      s += x;
      w->htot[ij] += x;
    }
    w->tot[k] = s;
    w->lntot[k] = (s > 0) ? LOG(s) : LOG0;
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
__inline static xdouble wham2_lnadd(xdouble a, xdouble b)
{
  xdouble c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + LOG(1 + EXP(-c));
}



/* save the density of states to file `fn` */
static int wham2_savelndos(wham2_t *w, const char *fn)
{
  hist2_t *hist = w->hist;
  FILE *fp;
  int i, j, ij, m = hist->m;
  xdouble x, y;

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
        fprintf(fp, "%g %g %g\n", (double) x, (double) y, (double) w->lndos[ij]);
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
  xdouble bx, by, x, y, xx, xy, yy, lnw, lnx, lny;
  xdouble lnz, slnx, slny, slnxx, slnxy, slnyy;
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
      lnx = LOG(x);
      for ( j = w->jmin; j < w->jmax; j++ ) {
        ij = i*m + j;
        if ( w->lndos[ij] <= LOG0 ) continue;
        y = (j + .5) * hist->dy;
        lny = LOG(y);

        lnw = w->lndos[ij] - bx * x - by * y;
        lnz = wham2_lnadd(lnz, lnw);
        slnx = wham2_lnadd(slnx, lnx + lnw);
        slny = wham2_lnadd(slny, lny + lnw);
        slnxx = wham2_lnadd(slnxx, lnx + lnx + lnw);
        slnxy = wham2_lnadd(slnxy, lnx + lny + lnw);
        slnyy = wham2_lnadd(slnyy, lny + lny + lnw);
      }
    }
    x = EXP(slnx - lnz);
    y = EXP(slny - lnz);
    xx = EXP(slnxx - lnz) - x * x;
    xy = EXP(slnxy - lnz) - x * y;
    yy = EXP(slnyy - lnz) - y * y;
    x += hist->xmin;
    y += hist->ymin;
    if ( fp != NULL )
      fprintf(fp, "%g %g %g %g %g %g %g %g\n",
          (double) bx, (double) by,
          (double) x, (double) y,
          (double) xx, (double) xy, (double) yy,
          (double) (lnz - bx * hist->xmin - by * hist->ymin));
  }
  if ( fp != NULL ) fclose(fp);
}



/* shift `lnz` such that `lnz[0] == 0` */
static void wham2_shift(xdouble *lnz, int nbeta)
{
  int i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



/* normalize by removing the center of mass motion */
static void wham2_rmcom(xdouble *arr,
    const xdouble *tot, int n)
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
static void wham2_normalize(xdouble *lnz, int nbeta, const void *ptr)
{
  const wham2_t *w = (const wham2_t *) ptr;

  if ( w->flags & WHAM2_RMCOM ) {
    wham2_rmcom(lnz, w->tot, nbeta);
  } else {
    wham2_shift(lnz, nbeta);
  }
}



/* estimate the partition function using the single histogram method */
static void wham2_estimatelnz(wham2_t *w, xdouble *lnz)
{
  hist2_t *hist = w->hist;
  int i, j, ij, n = hist->n, m = hist->m, nm = n * m;
  int nbeta = hist->rows, k, l, kk;
  xdouble dbx, dby, db, x, y, h, s, dlnz;
  xdouble dx = hist->dx, dy = hist->dy;

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
      x = hist->xmin + (i * 2 + 1) * dx / 2;
      for ( j = w->jmin; j < w->jmax; j++ ) {
        ij = i*m + j;

        h = hist->arr[k*nm + ij];
        if ( h <= 0 ) continue;
        y = hist->ymin + (j * 2 + 1) * dy / 2;

        s += h;
        dlnz = wham2_lnadd(dlnz, LOG(h) + dbx * x + dby * y);
      }
    }
    lnz[k] = lnz[kk] + (s > 0 ? LOG(s) - dlnz : 0);
  }

  wham2_normalize(lnz, nbeta, w);
}



/* compute the partition function from the density of states */
static void wham2_getlnz(wham2_t *w, xdouble *lnz)
{
  hist2_t *hist = w->hist;
  int i, j, ij, k, m = hist->m, nbeta = hist->rows;
  xdouble bx, by, x, y;
  xdouble dx = hist->dx, dy = hist->dy;

  for ( k = 0; k < nbeta; k++ ) {
    bx = w->bx[k];
    by = w->by[k];
    lnz[k] = LOG0;
    for ( i = w->imin; i < w->imax; i++ ) {
      x = hist->xmin + (i * 2 + 1) * dx / 2;
      for ( j = w->jmin; j < w->jmax; j++ ) {
        ij = i*m + j;
        if ( w->lndos[ij] <= LOG0) continue;
        y = hist->ymin + (j * 2 + 1) * dy / 2;

        /* res[k] is the new partition function */
        lnz[k] = wham2_lnadd(lnz[k], w->lndos[ij] - bx * x - by * y);
      }
    }
  }

  wham2_normalize(lnz, nbeta, w);
}



static xdouble wham2_step(wham2_t *w, xdouble *lnz, xdouble *res,
    xdouble damp)
{
  hist2_t *hist = w->hist;
  int i, j, ij, ijmin, n = hist->n, m = hist->m, nm = n * m;
  int k, nbeta = hist->rows;
  xdouble x, y, bx, by, lnden, err;
  xdouble dx = hist->dx, dy = hist->dy;

  ijmin = -1;
  for ( i = w->imin; i < w->imax; i++ ) {
    x = hist->xmin + (i * 2 + 1) * dx / 2;
    for ( j = w->jmin; j < w->jmax; j++ ) {
      y = hist->ymin + (j * 2 + 1) * dy / 2;
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
      w->lndos[ij] = LOG(w->lndos[ij]) - lnden;
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
    if ( FABS(res[i]) > err ) {
      err = FABS(res[i]);
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
static xdouble wham2_getlndos(wham2_t *w, xdouble *lnz,
    xdouble damp, int itmin, int itmax, xdouble tol, int verbose)
{
  int it;
  xdouble err, errp;
  clock_t t0, t1;

  t0 = clock();
  err = errp = 1e30;
  for ( it = 1; it <= itmax; it++ ) {
    err = wham2_step(w, lnz, w->res, damp);
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
  fprintf(stderr, "WHAM2 converged in %d steps, error %g, time %.4fs\n",
      it, (double) err, 1.0*(t1 - t0)/CLOCKS_PER_SEC);
  return err;
}



/* two-dimensional weighted histogram analysis method */
static xdouble wham2(hist2_t *hist,
    const xdouble *bx, const xdouble *by, xdouble *lnz,
    unsigned flags, xdouble damp,
    int itmin, int itmax, xdouble tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham2_t *w = wham2_open(bx, by, hist, flags);
  xdouble err;

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
#include "mdiis_xdbl.h"



static xdouble wham2_getres(void *w, xdouble *lnz, xdouble *res)
{
  return wham2_step((wham2_t *) w, lnz, res, 0);
}



static xdouble wham2_mdiis(hist2_t *hist,
    const xdouble *bx, const xdouble *by, xdouble *lnz,
    unsigned flags,
    int nbases, xdouble damp, int queue, xdouble threshold,
    int itmin, int itmax, xdouble tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham2_t *w = wham2_open(bx, by, hist, flags);
  xdouble err;

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



__inline static xdouble wham2x(hist2_t *hist,
    const xdouble *bx, const xdouble *by, xdouble *lnz,
    unsigned flags,
    xdouble damp, int nbases, int update_method, xdouble threshold,
    int itmin, int itmax, xdouble tol, int verbose,
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

