#ifndef WHAM_H__
#define WHAM_H__



/* Weighted histogram analysis method
 * this module requires `hist.h` */



#include <time.h>
#include "hist.h"



#ifndef LOG0
#define LOG0 -1e9
#endif



enum {
  WHAM_DIRECT = 0,
  WHAM_MDIIS = 1,
  WHAM_ST = 2,
  WHAM_UI = 3,
  WHAM_SCUI = 4,
  WHAM_NMETHODS
};

const char *wham_methods[] = {
  "Direct",
  "MDIIS",
  "ST",
  "UI",
  "SCUI",
  "WHAM_NMETHODS"
};



typedef struct {
  const double *beta; /* temperature array, reference */
  double *res; /* difference between new and old lnz */
  double *lndos; /* density of states */
  double *lntot; /* total number of visits to each temperature */
  const hist_t *hist; /* histograms, reference */
  double *sum, *ave, *var, *hnorm;
  double *ave1, *var1;
  int imin, imax;
} wham_t;



static wham_t *wham_open(const double *beta, const hist_t *hist)
{
  wham_t *w;
  int i, j, nbeta = hist->rows, n = hist->n;
  double h;

  xnew(w, 1);
  w->beta = beta;
  w->hist = hist;
  xnew(w->lntot, hist->rows);
  xnew(w->res, hist->rows);
  xnew(w->lndos, hist->n);
  xnew(w->sum, hist->rows);
  xnew(w->ave, hist->rows);
  xnew(w->var, hist->rows);
  xnew(w->hnorm, hist->rows);
  xnew(w->ave1, hist->rows);
  xnew(w->var1, hist->rows);

  /* compute the total */
  for ( j = 0; j < nbeta; j++ ) {
    double x = 0;
    for ( i = 0; i < n; i++ )
      x += hist->arr[j*n+i];
    w->lntot[j] = (x > 0) ? log(x) : LOG0;
  }

  /* determine the boundaries */
  /* find imin */
  for ( i = 0; i < n; i++ ) {
    for ( h = 0, j = 0; j < nbeta; j++ ) {
      h += hist->arr[j * n + i];
    }
    if ( h > 0 ) break;
  }
  w->imin = i;

  /* find imax */
  for ( i = n - 1; i >= 0; i-- ) {
    for ( h = 0, j = 0; j < nbeta; j++ ) {
      h += hist->arr[j * n + i];
    }
    if ( h > 0 ) break;
  }
  w->imax = i + 1;

  fprintf(stderr, "i: [%d, %d)\n", w->imin, w->imax);
  return w;
}



static void wham_close(wham_t *w)
{
  free(w->res);
  free(w->lndos);
  free(w->lntot);
  free(w->sum);
  free(w->ave);
  free(w->var);
  free(w->hnorm);
  free(w->ave1);
  free(w->var1);
  free(w);
}



/* log(exp(a) + exp(b)) */
__inline static double wham_lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a - b) > 50.0) ? a : a + log(1 + exp(-c));
}



/* save the density of states to file `fn` */
static int wham_savelndos(wham_t *w, const char *fn)
{
  FILE *fp;
  int i;
  double emin = w->hist->xmin, de = w->hist->dx;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  for ( i = w->imin; i < w->imax; i++ )
    if ( w->lndos[i] > LOG0 )
      fprintf(fp, "%g %g\n", emin + (i+.5)*de, w->lndos[i]);
  fclose(fp);
  return 0;
}



/* compute the energy and heat capacity from the WHAM */
static void wham_getav(wham_t *w, const char *fn)
{
  const hist_t *hist = w->hist;
  double T, b, e, ee, lne, lnz, slne, slnee, lnw;
  double de = hist->dx, emin = hist->xmin;
  int i, j, nbeta = hist->rows;
  FILE *fp = NULL;

  if ( fn != NULL && (fp = fopen(fn, "w")) == NULL )
    fprintf(stderr, "cannot write %s\n", fn);

  for ( j = 0; j < nbeta; j++ ) {
    b = w->beta[j];
    T = 1./b;
    lnz = slne = slnee = LOG0;
    for ( i = w->imin; i < w->imax; i++ ) {
      if ( w->lndos[i] <= LOG0 ) continue;
      /* note: we do not add emin here for it may lead to
       * a negative energy whose logarithm is undefined */
      e = (i + .5) * de;
      lne = log(e);
      lnw = w->lndos[i] - b * e;
      lnz = wham_lnadd(lnz, lnw);
      slne = wham_lnadd(slne, lne + lnw);
      slnee = wham_lnadd(slnee, 2*lne + lnw);
    }
    e = exp(slne - lnz);
    ee = exp(slnee - lnz) - e * e;
    e += emin;
    if ( fp != NULL )
      fprintf(fp, "%g %g %g %g\n", T, e, ee/(T*T), lnz - b * emin);
  }
  if ( fp != NULL ) fclose(fp);
}



/* estimate the partition function using the single histogram method */
static void wham_estimatelnz(wham_t *w, double *lnz)
{
  const hist_t *hist = w->hist;
  int i, j, n = hist->n, nbeta = hist->rows;
  double db, e, h, s, dlnz;

  lnz[0] = 0;
  for ( j = 1; j < nbeta; j++ ) {
    /* estimate the free energy different between
     * the pair j - 1 and j */
    db = w->beta[j] - w->beta[j - 1];
    s = 0;
    dlnz = LOG0;
    /* Z(j-1) / Z(j) = < exp( [beta(j) - beta(j-1)] E ) >_j
     *    Sum_E h_j(E) exp( [beta(j) - beta(j-1)] E )
     * = ---------------------------------------------
     *               Sum_E h_j(E)
     **/
    for ( i = w->imin; i < w->imax; i++ ) {
      h = hist->arr[j*n + i];
      if ( h <= 0 ) continue;
      e = hist->xmin + (i + .5) * hist->dx;
      s += h;
      dlnz = wham_lnadd(dlnz, log(h) + db * e);
    }
    lnz[j] = lnz[j - 1] + (s > 0 ? log(s) - dlnz : 0);
  }
}



static void wham_normalize(double *lnz, int nbeta)
{
  int i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



/* compute the partition function from the density of states */
static void wham_getlnz(wham_t *w, double *lnz)
{
  const hist_t *hist = w->hist;
  int i, j, nbeta = hist->rows;
  double e;

  for ( j = 0; j < nbeta; j++ ) {
    lnz[j] = LOG0;
    for ( i = w->imin; i < w->imax; i++ ) {
      if ( w->lndos[i] <= LOG0 ) continue;
      e = hist->xmin + (i + .5) * hist->dx;
      lnz[j] = wham_lnadd(lnz[j], w->lndos[i] - w->beta[j] * e);
    }
  }
  wham_normalize(lnz, nbeta);
}



static double wham_step(wham_t *w, double *lnz, double *res,
    double damp)
{
  const hist_t *hist = w->hist;
  int i, j, imin, n = hist->n, nbeta = hist->rows;
  double x, num, lnden, e, emin = hist->xmin, de = hist->dx, err;

  imin = -1;
  for ( i = w->imin; i < w->imax; i++ ) {
    num = 0;
    lnden = LOG0;
    e = emin + (i + .5) * de;
    /*        num           Sum_j h_j(i)
     * dos = ----- = ------------------------------------------
     *        den     Sum_j tot_j exp(-beta_j * e) / Z_j
     * */
    for ( j = 0; j < nbeta; j++ ) {
      x = hist->arr[j*n + i];
      //if ( x <= 0 ) continue;
      num += x;
      lnden = wham_lnadd(lnden, w->lntot[j] - w->beta[j] * e - lnz[j]);
    }
    if ( num > 0 ) {
      w->lndos[i] = log(num) - lnden;
      if ( imin < 0 ) imin = i;
    } else {
      w->lndos[i] = LOG0;
    }
  }

  /* shift the baseline of the density of states */
  for ( x = w->lndos[imin], i = 0; i < n; i++ )
    if ( w->lndos[i] > LOG0 )
      w->lndos[i] -= x;

  /* refresh the partition function */
  wham_getlnz(w, res);

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
    wham_normalize(lnz, nbeta);
  }

  return err;
}



/* iteratively compute the logarithm of the density of states
 * using the weighted histogram method */
static double wham_getlndos(wham_t *w, double *lnz,
    double damp, int itmin, int itmax, double tol, int verbose)
{
  int it;
  double err, errp;
  clock_t t0, t1;

  t0 = clock();
  err = errp = 1e30;
  for ( it = 0; it < itmax; it++ ) {
    err = wham_step(w, lnz, w->res, damp);
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
  fprintf(stderr, "WHAM converged in %d steps, error %g, time %.4fs\n",
      it, err, 1.0*(t1 - t0)/CLOCKS_PER_SEC);
  return err;
}



/* weighted histogram analysis method */
static double wham(const hist_t *hist, const double *beta, double *lnz,
    double damp, int itmin, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);
  double err;

  wham_estimatelnz(w, lnz);
  err = wham_getlndos(w, lnz,
      damp, itmin, itmax, tol, verbose);
  if ( fnlndos ) {
    wham_savelndos(w, fnlndos);
  }
  wham_getav(w, fneav);
  wham_close(w);
  return err;
}



/* non-iteratively compute the logarithm of the density of states
 * using the statistical-temperature WHAM
 * if `gauss`, use Gaussian correlation for missing data */
static void stwham_getlndos(wham_t *w)
{
  const hist_t *hist = w->hist;
  int i, j, imin, imax, n = hist->n, nbeta = hist->rows;
  double x, y, y2, tot, hn, hp, stbeta, de = hist->dx;
  clock_t t0, t1;

  t0 = clock();

  imin = w->imin;
  imax = w->imax;
  for ( i = imin; i < imax; i++ ) {
    tot = 0;
    stbeta = 0;
    for ( j = 0; j < nbeta; j++ ) {
      x = hist->arr[j*n + i];
      tot += x;
      /* compute the statistical temperature from copy j
       *            sum_j [(n_j)'(E) + n_j(E) beta_j]
       * beta(E) = -----------------------------------
       *                     sum_k n_k(E)
       * y is the numerator */
      y = x * w->beta[j];
      y2 = 0;
      if ( i - 1 >= imin && i + 1 < imax ) {
        hn = hist->arr[j*n + i + 1];
        hp = hist->arr[j*n + i - 1];
        /* compute (n_j)'(E) + n_j(E) beta */
        y2 = (hn - hp) / (2 * de);
      }
      stbeta += y + y2;
    }

    /* lndos currently holds the statistical temperature */
    if ( tot > 0 ) {
      w->lndos[i] = stbeta / tot;
    } else {
      w->lndos[i] = 0;
      continue;
    }
  }

  /* integrate the statistical temperature
   * to get the density of states */
  x = y = 0;
  for ( i = imin; i < imax; i++ ) {
    stbeta = w->lndos[i];
    w->lndos[i] = y + 0.5 * (x + stbeta) * de;
    x = stbeta; /* previous temperature */
    y = w->lndos[i]; /* previous lndos */
  }

  for ( i = 0; i < imin; i++ ) {
    w->lndos[i] = LOG0;
  }
  for ( i = w->imax; i < n; i++ ) {
    w->lndos[i] = LOG0;
  }

  t1 = clock();
  fprintf(stderr, "ST-WHAM completed, time %.4fs\n",
      1.0*(t1 - t0)/CLOCKS_PER_SEC);
}



/* statistical temperature weighted histogram analysis method */
static double stwham(const hist_t *hist, const double *beta, double *lnz,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);

  stwham_getlndos(w);
  wham_getlnz(w, lnz);
  if ( fnlndos ) {
    wham_savelndos(w, fnlndos);
  }
  wham_getav(w, fneav);
  wham_close(w);
  return 0.0;
}



/* compute histogram averages */
static void wham_gethave(wham_t *w)
{
  const hist_t *hist = w->hist;
  int i, j, n = hist->n, nbeta = hist->rows;
  double tot, sx, sxx, x, y, de = hist->dx;

  /* compute averages and variance */
  for ( j = 0; j < nbeta; j++ ) {
    tot = sx = sxx = 0;
    for ( i = w->imin; i < w->imax; i++ ) {
      y = hist->arr[j*n + i];
      x = hist->xmin + ( i + 0.5 ) * de;
      tot += y;
      sx += x * y;
      sxx += x * x * y;
    }
    sx /= tot;
    sxx = sxx / tot - sx * sx;
    if ( tot > 1 ) {
      sxx *= tot / (tot - 1);
    }
    w->sum[j] = tot;
    w->ave[j] = sx;
    w->var[j] = sxx;
  }
}



/* non-iteratively compute the logarithm of the density of states
 * using umbrella integration */
static void umbint_getlndos(wham_t *w)
{
  const hist_t *hist = w->hist;
  int i, j, imin, imax, n = hist->n, nbeta = hist->rows;
  double x, y, tot, stbeta, de = hist->dx;
  clock_t t0, t1;

  t0 = clock();

  /* compute averages and variance */
  wham_gethave(w);
  for ( j = 0; j < nbeta; j++ ) {
    w->hnorm[j] = w->sum[j] * de / sqrt(2 * M_PI * w->var[j]);
  }

  imin = w->imin;
  imax = w->imax;
  for ( i = imin; i < imax; i++ ) {
    double ei = hist->xmin + (i + 0.5) * de;

    tot = 0;
    stbeta = 0;
    for ( j = 0; j < nbeta; j++ ) {
      /* use Gaussian approximation */
      double De = w->ave[j] - ei;
      x = w->hnorm[j] * exp( -0.5 * De * De / w->var[j] );
      tot += x;
      stbeta += x * (w->beta[j] + De / w->var[j]);
    }

    /* lndos currently holds the statistical temperature */
    w->lndos[i] = (tot > 0) ? stbeta / tot : 0;
  }

  /* integrate the statistical temperature
   * to get the density of states */
  x = y = 0;
  for ( i = imin; i < imax; i++ ) {
    stbeta = w->lndos[i];
    w->lndos[i] = y + 0.5 * (x + stbeta) * de;
    x = stbeta; /* previous temperature */
    y = w->lndos[i]; /* previous lndos */
  }

  for ( i = 0; i < imin; i++ ) {
    w->lndos[i] = LOG0;
  }
  for ( i = w->imax; i < n; i++ ) {
    w->lndos[i] = LOG0;
  }

  t1 = clock();
  fprintf(stderr, "ST-WHAM completed, time %.4fs\n",
      1.0*(t1 - t0)/CLOCKS_PER_SEC);
}



/* umbrella integration */
static double umbint(const hist_t *hist, const double *beta, double *lnz,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);

  umbint_getlndos(w);
  wham_getlnz(w, lnz);
  if ( fnlndos ) {
    wham_savelndos(w, fnlndos);
  }
  wham_getav(w, fneav);
  wham_close(w);
  return 0.0;
}



/* iteratively compute the logarithm of the density of states
 * using self-consistent umbrella integration */
static void scui_getlndos(wham_t *w,
    double damp, int itmin, int itmax, double tol, int verbose)
{
  const hist_t *hist = w->hist;
  int it, i, j, imin, imax, n = hist->n, nbeta = hist->rows;
  double x, y, tot, stbeta, de = hist->dx, err, errp;
  clock_t t0, t1;

  t0 = clock();
  err = errp = 1e30;
  /* compute averages and variance */
  wham_gethave(w);
  for ( j = 0; j < nbeta; j++ ) {
    w->ave1[j] = w->ave[j];
    w->var1[j] = w->var[j];
  }

  imin = w->imin;
  imax = w->imax;
  for ( it = 0; it < itmax; it++ ) {
    for ( j = 0; j < nbeta; j++ ) {
      w->hnorm[j] = w->sum[j] * de / sqrt(2 * M_PI * w->var1[j]);
    }

    err = 0;
    for ( i = imin; i < imax; i++ ) {
      double ei = hist->xmin + (i + 0.5) * de;

      tot = 0;
      stbeta = 0;
      for ( j = 0; j < nbeta; j++ ) {
        /* use Gaussian approximation */
        double De = w->ave1[j] - ei;
        x = w->hnorm[j] * exp( -0.5 * De * De / w->var1[j] );
        tot += x;
        stbeta += x * (w->beta[j] + De / w->var1[j]);
      }

      /* lndos currently holds the statistical temperature */
      w->lndos[i] = (tot > 0) ? stbeta / tot : 0;
    }

    /* integrate the statistical temperature
     * to get the density of states */
    x = y = 0;
    for ( i = imin; i < imax; i++ ) {
      stbeta = w->lndos[i];
      w->lndos[i] = y + 0.5 * (x + stbeta) * de;
      x = stbeta; /* previous temperature */
      y = w->lndos[i]; /* previous lndos */
    }

    int jj = -1;
    /* update the averages and variances */
    for ( j = 0; j < nbeta; j++ ) {
      double lnz, slne, slnee, e, lne, lnw, ee;
      lnz = slne = slnee = LOG0;
      for ( i = imin; i < imax; i++ ) {
        /* note: we do not add emin here for it may lead to
         * a negative energy whose logarithm is undefined */
        e = (i + .5) * de;
        lne = log(e);
        lnw = w->lndos[i] - w->beta[j] * e;
        lnz = wham_lnadd(lnz, lnw);
        slne = wham_lnadd(slne, lne + lnw);
        slnee = wham_lnadd(slnee, 2*lne + lnw);
      }
      e = exp(slne - lnz);
      ee = exp(slnee - lnz) - e * e;
      e += hist->xmin;
      if ( verbose >= 2 ) {
        fprintf(stderr, "j %d: ave %g vs %g (%g), var %g vs %g (%g)\n", j, e, w->ave1[j], w->ave[j], ee, w->var1[j], w->var[j]);
      }

      x = w->ave[j] - e;
      if ( fabs(x) > err ) { jj = j; err = fabs(x); }
      w->ave1[j] += damp * x;

      x = log( w->var[j] / ee );
      if ( fabs(x) > err ) { jj = j; err = fabs(x); }
      w->var1[j] *= exp(damp * x);
    }

    if ( verbose ) {
      fprintf(stderr, "it %d, err %g -> %g, jj %d\n",
          it, errp, err, jj);
      if ( verbose >= 2 ) getchar();
    }
    if ( err < tol && it > itmin ) {
      break;
    }
    if ( err > errp ) break;
    errp = err;
  }

  for ( i = 0; i < imin; i++ ) {
    w->lndos[i] = LOG0;
  }
  for ( i = w->imax; i < n; i++ ) {
    w->lndos[i] = LOG0;
  }

  t1 = clock();
  fprintf(stderr, "ST-WHAM completed, time %.4fs\n",
      1.0*(t1 - t0)/CLOCKS_PER_SEC);
}



/* self-consistent umbrella integration */
static double scui(const hist_t *hist, const double *beta, double *lnz,
    double damp, int itmin, int itmax, double tol,
    int verbose, const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);

  scui_getlndos(w, damp, itmin, itmax, tol, verbose);
  wham_getlnz(w, lnz);
  if ( fnlndos ) {
    wham_savelndos(w, fnlndos);
  }
  wham_getav(w, fneav);
  wham_close(w);
  return 0.0;
}



#ifdef ENABLE_MDIIS
/* MDIIS method */
#include "mdiis.h"



static double wham_getres(void *w, double *lnz, double *res)
{
  return wham_step((wham_t *) w, lnz, res, 0);
}



static double wham_mdiis(const hist_t *hist, const double *beta, double *lnz,
    int nbases, double damp, int queue, double threshold,
    int itmin, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);
  double err;

  wham_estimatelnz(w, lnz);
  err = iter_mdiis(lnz, hist->rows,
      wham_getres, wham_normalize, w,
      nbases, damp, queue, threshold,
      itmin, itmax, tol, verbose);
  if ( fnlndos ) wham_savelndos(w, fnlndos);
  wham_getav(w, fneav);
  wham_close(w);
  return err;
}



#endif /* ENABLE_MDIIS */



static double whamx(const hist_t *hist, const double *beta, double *lnz,
    double damp, int nbases, int update_method, double threshold,
    int itmin, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav, int method)
{
  if ( method == WHAM_DIRECT ) {
    return wham(hist, beta, lnz,
        damp, itmin, itmax, tol,
        verbose, fnlndos, fneav);
  } else if ( method == WHAM_ST ) {
    return stwham(hist, beta, lnz, fnlndos, fneav);
  } else if ( method == WHAM_UI ) {
    return umbint(hist, beta, lnz, fnlndos, fneav);
  } else if ( method == WHAM_SCUI ) {
    return scui(hist, beta, lnz,
        damp, itmin, itmax, tol,
        verbose, fnlndos, fneav);
#ifdef ENABLE_MDIIS
  } else if ( method == WHAM_MDIIS ) {
    return wham_mdiis(hist, beta, lnz,
        nbases, damp, update_method, threshold,
        itmin, itmax, tol, verbose,
        fnlndos, fneav);
#endif
  }

  return 0;
}



#endif /* WHAM_H__ */

