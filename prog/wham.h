#ifndef WHAM_H__
#define WHAM_H__



/* Weighted histogram analysis method
 * this module requires `hist.h` */



#include "hist.h"



#ifndef LOG0
#define LOG0 -1e9
#endif



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

  /* refresh the partition function */
  for ( j = 0; j < nbeta; j++ ) {
    for ( res[j] = LOG0, i = 0; i < n; i++ ) {
      if ( w->lndos[i] <= LOG0) continue;
      e = hist->xmin + (i + .5) * hist->dx;
      res[j] = wham_lnadd(res[j], w->lndos[i] - w->beta[j] * e);
    }
  }
  for ( x = res[0], j = 0; j < nbeta; j++ )
    res[j] -= x; /* shift the baseline */

  for ( err = 0, j = 0; j < nbeta; j++ ) {
    res[j] -= lnz[j];
    if ( fabs(res[j]) > err ) err = fabs(res[j]);
    if ( update ) lnz[j] += res[j];
  }

  return err;
}



/* iteratively compute the logarithm of the density of states
 * using the weighted histogram method */
static double wham_getlndos(wham_t *w, double *lnz,
    int itmax, double tol, int verbose)
{
  int it;
  double err, errp = 1e30;

  for ( it = 1; it <= itmax; it++ ) {
    err = wham_step(w, lnz, w->res, 1);
    if ( verbose ) {
      fprintf(stderr, "it %d, err %g -> %g\n",
          it, errp, err);
    }
    if ( err < tol ) {
      break;
    }
    errp = err;
  }

  fprintf(stderr, "WHAM converged at step %d, error %g\n", it, err);
  return err;
}



/* weighted histogram analysis method */
static double wham(hist_t *hist, const double *beta, double *lnz,
    int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);
  double err;

  wham_estimatelnz(w, lnz);
  err = wham_getlndos(w, lnz, itmax, tol, verbose);
  if ( fnlndos ) {
    wham_savelndos(w, fnlndos);
  }
  wham_getav(w, fneav);
  wham_close(w);
  return err;
}



#ifdef WHAM_MDIIS
/* MDIIS method */
#include "mdiis.h"



static double wham_getres(void *w, double *lnz, double *res)
{
  return wham_step((wham_t *) w, lnz, res, 0);
}



static double wham_mdiis(hist_t *hist, const double *beta, double *lnz,
    int nbases, double damp, int itmax, double tol, int verbose,
    const char *fnlndos, const char *fneav)
{
  wham_t *w = wham_open(beta, hist);
  double err;

  wham_estimatelnz(w, lnz);
  err = iter_mdiis(lnz, hist->rows, wham_getres, w,
      nbases, damp, itmax, tol, verbose);
  if ( fnlndos ) wham_savelndos(w, fnlndos);
  wham_getav(w, fneav);
  wham_close(w);
  return err;
}

#endif /* WHAM_MDIIS */



#endif /* WHAM_H__ */

