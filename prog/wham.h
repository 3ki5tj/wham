#ifndef WHAM_H__
#define WHAM_H__



/* Weighted histogram analysis method
 * this module requires `hist.h` */



#include <time.h>
#include "hist.h"



#ifndef LOG0
#define LOG0 -1e9
#endif



enum { WHAM_DIRECT = 0, WHAM_MDIIS = 1, WHAM_ST = 2, WHAM_NMETHODS };
const char *wham_methods[] = {"Direct", "MDIIS", "ST", "WHAM_NMETHODS"};



typedef struct {
  const double *beta; /* temperature array, reference */
  double *res; /* difference between new and old lnz */
  double *lndos; /* density of states */
  double *lntot; /* total number of visits to each temperature */
  hist_t *hist; /* histograms, reference */
} wham_t;



static wham_t *wham_open(const double *beta, hist_t *hist)
{
  wham_t *w;
  int i, j, nbeta = hist->rows, n = hist->n;

  xnew(w, 1);
  w->beta = beta;
  w->hist = hist;
  xnew(w->lntot, hist->rows);
  xnew(w->res, hist->rows);
  xnew(w->lndos, hist->n);

  /* compute the total */
  for ( j = 0; j < nbeta; j++ ) {
    double x = 0;
    for ( i = 0; i < n; i++ )
      x += hist->arr[j*n+i];
    w->lntot[j] = (x > 0) ? log(x) : LOG0;
  }

  return w;
}



static void wham_close(wham_t *w)
{
  free(w->res);
  free(w->lndos);
  free(w->lntot);
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
  int i, n = w->hist->n;
  double emin = w->hist->xmin, de = w->hist->dx;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  for ( i = 0; i < n; i++ )
    if ( w->lndos[i] > LOG0 )
      fprintf(fp, "%g %g\n", emin + (i+.5)*de, w->lndos[i]);
  fclose(fp);
  return 0;
}



/* compute the energy and heat capacity from the WHAM */
static void wham_getav(wham_t *w, const char *fn)
{
  hist_t *hist = w->hist;
  double T, b, e, ee, lne, lnz, slne, slnee, lnw;
  double de = hist->dx, emin = hist->xmin;
  int i, j, n = hist->n, nbeta = hist->rows;
  FILE *fp = NULL;

  if ( fn != NULL && (fp = fopen(fn, "w")) == NULL )
    fprintf(stderr, "cannot write %s\n", fn);

  for ( j = 0; j < nbeta; j++ ) {
    b = w->beta[j];
    T = 1./b;
    lnz = slne = slnee = LOG0;
    for ( i = 0; i < n; i++ ) {
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
  hist_t *hist = w->hist;
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
    for ( i = 0; i < n; i++ ) {
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
  hist_t *hist = w->hist;
  int i, j, n = hist->n, nbeta = hist->rows;
  double e;

  for ( j = 0; j < nbeta; j++ ) {
    for ( lnz[j] = LOG0, i = 0; i < n; i++ ) {
      if ( w->lndos[i] <= LOG0) continue;
      e = hist->xmin + (i + .5) * hist->dx;
      lnz[j] = wham_lnadd(lnz[j], w->lndos[i] - w->beta[j] * e);
    }
  }
  wham_normalize(lnz, nbeta);
}



static double wham_step(wham_t *w, double *lnz, double *res, int update)
{
  hist_t *hist = w->hist;
  int i, j, imin, n = hist->n, nbeta = hist->rows;
  double x, num, lnden, e, emin = hist->xmin, de = hist->dx, err;

  imin = -1;
  for ( i = 0; i < n; i++ ) {
    num = 0;
    lnden = LOG0;
    e = emin + (i + .5) * de;
    /*        num           Sum_j h_j(i)
     * dos = ----- = ------------------------------------------
     *        den     Sum_j tot_j exp(-beta_j * e) / Z_j
     * */
    for ( j = 0; j < nbeta; j++ ) {
      x = hist->arr[j*n + i];
      if ( x <= 0 ) continue;
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

  wham_getlnz(w, res);
  for ( err = 0, j = 0; j < nbeta; j++ ) {
    x = res[j];
    res[j] = x - lnz[j];
    if ( fabs(res[j]) > err ) err = fabs(res[j]);
    if ( update ) lnz[j] = x;
  }

  return err;
}



/* iteratively compute the logarithm of the density of states
 * using the weighted histogram method */
static double wham_getlndos(wham_t *w, double *lnz,
    int itmax, double tol, int itmin, int verbose)
{
  int it;
  double err, errp;
  clock_t t0, t1;

  t0 = clock();
  err = errp = 1e30;
  for ( it = 0; it < itmax; it++ ) {
    err = wham_step(w, lnz, w->res, 1);
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
static double wham(hist_t *hist, const double *beta, double *lnz,
    int itmax, double tol, int itmin, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);
  double err;

  wham_estimatelnz(w, lnz);
  err = wham_getlndos(w, lnz, itmax, tol, itmin, verbose);
  if ( fnlndos ) {
    wham_savelndos(w, fnlndos);
  }
  wham_getav(w, fneav);
  wham_close(w);
  return err;
}



/* non-iteratively compute the logarithm of the density of states
 * using the statistical-temperature WHAM */
static void stwham_getlndos(wham_t *w)
{
  hist_t *hist = w->hist;
  int i, j, imin, imax, n = hist->n, nbeta = hist->rows;
  double x, y, tot, hn, hp, stbeta, de = hist->dx;
  clock_t t0, t1;

  t0 = clock();

  imin = imax = -1;
  for ( i = 0; i < n; i++ ) {
    stbeta = 0;
    tot = 0;
    for ( j = 0; j < nbeta; j++ ) {
      x = hist->arr[j*n + i];
      if ( x <= 0 ) continue;
      tot += x;
      /* compute the statistical temperature from copy j
       *            sum_j [(n_j)'(E) + n_j(E) beta_j]
       * beta(E) = -----------------------------------
       *                     sum_k n_k(E)
       * y is the denominator */
      y = x * w->beta[j];
      if ( i - 1 >= 0 && i + 1 < n ) {
        hn = hist->arr[j*n + i + 1];
        hp = hist->arr[j*n + i - 1];
        /* compute (n_j)'(E) + n_j(E) beta */
        y += (hn - hp) / (2 * de);
      }
      stbeta += y;
    }

    /* lndos currently holds the statistical temperature */
    if ( tot > 0 ) {
      w->lndos[i] = stbeta / tot;
    } else {
      w->lndos[i] = 0;
      continue;
    }
  
    if ( imin < 0 ) imin = i;
    if ( i > imax ) imax = i;
  }

  /* integrate the statistical temperature
   * to get the density of states */
  x = y = 0;
  for ( i = imin; i <= imax; i++ ) {
    stbeta = w->lndos[i];
    w->lndos[i] = y + 0.5 * (x + stbeta) * de;
    x = stbeta; /* previous temperature */
    y = w->lndos[i]; /* previous lndos */
  }

  for ( i = 0; i < imin; i++ ) {
    w->lndos[i] = LOG0;
  }
  for ( i = imax + 1; i < n; i++ ) {
    w->lndos[i] = LOG0;
  }

  t1 = clock();
  fprintf(stderr, "ST-WHAM completed, time %.4fs\n",
      1.0*(t1 - t0)/CLOCKS_PER_SEC);
}



/* statistical temperature weighted histogram analysis method */
static double stwham(hist_t *hist, const double *beta, double *lnz,
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



#ifdef ENABLE_MDIIS
/* MDIIS method */
#include "mdiis.h"



static double wham_getres(void *w, double *lnz, double *res)
{
  return wham_step((wham_t *) w, lnz, res, 0);
}



static double wham_mdiis(hist_t *hist, const double *beta, double *lnz,
    int nbases, double damp, int queue, double threshold,
    int itmax, double tol, int itmin, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);
  double err;

  wham_estimatelnz(w, lnz);
  err = iter_mdiis(lnz, hist->rows,
      wham_getres, wham_normalize, w,
      nbases, damp, queue, threshold,
      itmax, tol, itmin, verbose);
  if ( fnlndos ) wham_savelndos(w, fnlndos);
  wham_getav(w, fneav);
  wham_close(w);
  return err;
}



#endif /* ENABLE_MDIIS */



static double whamx(hist_t *hist, const double *beta, double *lnz,
    int nbases, double damp, int update_method, double threshold,
    int itmax, double tol, int itmin, int verbose,
    const char *fnlndos, const char *fneav, int method)
{
  if ( method == WHAM_DIRECT ) {
    return wham(hist, beta, lnz, itmax, tol, itmin, verbose,
        fnlndos, fneav);
  } else if ( method == WHAM_ST ) {
    return stwham(hist, beta, lnz, fnlndos, fneav);
#ifdef ENABLE_MDIIS
  } else if ( method == WHAM_MDIIS ) {
    return wham_mdiis(hist, beta, lnz,
        nbases, damp, update_method, threshold,
        itmax, tol, itmin, verbose,
        fnlndos, fneav);
#endif
  }

  return 0;
}



#endif /* WHAM_H__ */

