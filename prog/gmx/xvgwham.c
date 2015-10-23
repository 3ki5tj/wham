/* WHAM for a set of xvg files */
#define ENABLE_MDIIS
#include "../wham_xdbl.h"
#include <time.h>
#include "../mtrand.h"
#define MTRAND
#include "xvg.h"
#define WHAM
#include "../whammodel.h"
#include "lsutil.h"





/* save autocorrelation times */
static int saveact(xvg_t **xvg, int nbet,
    double actmax, double acmin, const char *fnact)
{
  int i, j;
  xvg_t *xvg0;
  double act[10];
  FILE *fp;

  if ( fnact == NULL ) {
    return -1;
  } else if ( (fp = fopen(fnact, "w")) == NULL ) {
    fprintf(stderr, "cannot write file %s\n", fnact);
    return -1;
  }

  fprintf(fp, "# %d %g\n", nbet, xvg[0]->dx);

  for ( i = 0; i < nbet; i++ ) {
    xvg0 = xvg[i];
    xvg_act(xvg0, act, actmax, acmin, NULL);
    for ( j = 0; j < xvg0->m; j++ ) {
      fprintf(fp, " %g", act[j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}



/* load autocorrelation times in terms of number of frames */
static int loadact(const char *fnact, int nbet, double *tcorr)
{
  int i, err = 0;
  double dx;
  char s[1024];
  FILE *fp;

  if ( fnact == NULL ) {
    return -1;
  } else if ( (fp = fopen(fnact, "r")) == NULL ) {
    fprintf(stderr, "cannot load file %s\n", fnact);
    return -1;
  }

  for ( i = 0; i < nbet; i++ ) {
    tcorr[i] = 0;
  }

  if ( fgets(s, sizeof s, fp) == NULL ) {
    err = -1;
    goto END;
  }
  sscanf(s + 1, "%d%lf", &i, &dx);

  for ( i = 0; i < nbet; i++ ) {
    if ( fgets(s, sizeof s, fp) == NULL ) {
      fprintf(stderr, "cannot load entry %d from %s\n", i, fnact);
      err = -1;
      break;
    }
    sscanf(s, "%lf", &tcorr[i]);
    tcorr[i] /= dx;
    //fprintf(stderr, "tcorr %d: %g %g\n", i, tcorr[i], tcorr[i] * dx);
  }

END:
  fclose(fp);
  return err;
}



/* construct the histogram
 * optionally compute autocorrelation times
 * and save to `fnact` */
static hist_t *mkhist(const char *fnls,
    xdouble **beta, double de,
    const char *fnhis, double radd,
    double actmax, double acmin, const char *fnact)
{
  hist_t *hs;
  int i, j, nbeta;
  char **fns;
  xvg_t **xvg = NULL;
  double xmin[2], xmax[2], emin = 1e30, emax = -1e30;

  /* scramble the random number seed */
  mtscramble( time(NULL) + clock() );

  if ( (fns = getls(fnls, &nbeta, beta)) == NULL ) {
    return NULL;
  }

  xnew(xvg, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    if ( (xvg[i] = xvg_load(fns[i], radd)) == NULL ) {
      exit(1);
    }
    xvg_minmax(xvg[i], xmin, xmax);
    if ( xmin[0] < emin ) {
      emin = xmin[0];
    }
    if ( xmax[0] > emax ) {
      emax = xmax[0];
    }
  }

  emin = ((int) (emin / de) - 1) * de;
  emax = ((int) (emax / de) + 1) * de;

  hs = hist_open(nbeta, emin, emax, de);

  for ( i = 0; i < nbeta; i++ ) {
    for ( j = 0; j < xvg[i]->n; j++ ) {
      hist_add1(hs, i, xvg[i]->y[0][j], 1.0, HIST_VERBOSE);
    }
  }

  hist_save(hs, fnhis, HIST_ADDAHALF | HIST_NOZEROES | HIST_VERBOSE);

  /* save autocorrelation time */
  saveact(xvg, nbeta, actmax, acmin, fnact);

  for ( i = 0; i < nbeta; i++ ) {
    xvg_close(xvg[i]);
    free(fns[i]);
  }
  free(xvg);
  free(fns);

  return hs;
}



/* bootstrapping individual histograms
 * The autocorrelation time is `tau`
 * which is measured in terms of the number of frames */
static hist_t *hist_bootstrap(hist_t *hs0, double *tau)
{
  hist_t *hs;
  int i, imin, imax, k, r, n, il, ir;
  double x, tot;
  double *arr, *cnt;
  double gam;

  hs = hist_open(hs0->rows, hs0->xmin,
      hs0->xmin + hs0->n * hs0->dx, hs0->dx);

  n = hs->n;
  xnew(cnt, n + 1);

  for ( r = 0; r < hs->rows; r++ ) {
    /* bootstrap for histogram r */
    arr = hs0->arr + r * hs->n;

    /* count the total for histogram r */
    tot = 0;
    imin = n;
    imax = 0;
    cnt[0] = 0;
    for ( i = 0; i < n; i++ ) {
      x = arr[i];
      cnt[i + 1] = cnt[i] + x;
      if ( x > 0 ) {
        if ( i < imin ) imin = i;
        if ( i + 1 > imax ) imax = i + 1;
      }
    }
    tot = cnt[n];

    /* bootstrapping based on Hub 2010 */
    x = randgaus();
    gam = ( tau != NULL && tau[r] > 0 ) ? exp(-1/tau[r]) : 0;

    for ( i = 0, k = 0; k < tot; k++ ) {
      /* randomly change i */
      if ( k == 0 || rand01() >= gam ) {
        /* translate y to the histogram frame */
        x = tot * rand01();
        /* find the bin i containing x using binary search
         * that is cnt[i] < x < cnt[i+1] */
        il = imin;
        ir = imax;
        while ( il < ir - 1 ) {
          i = (il + ir + 1) / 2;
          if ( x > cnt[i] ) {
            il = i;
          } else {
            ir = i;
          }
        }
        i = il;
      }

      hs->arr[r * n + i] += 1;
    }
  }

  free(cnt);
  return hs;
}



/* bootstrapping individual histograms
 * This algorithm is based on
 *
 * g_wham -- A free weighted histogram analysis
 * implementation including robust error and
 * autocorrelation estimates
 * Jochen S. Hub, Bert L. de Groot, and David van der Spoel
 * J. Chem. Theory Comput. 2010 6, 3713-3720
 *
 * The autocorrelation time is `tau`
 * which is measured in terms of the number of frames */
__inline static hist_t *hist_bootstrap_hub(hist_t *hs0, double *tau)
{
  hist_t *hs;
  int i, imin, imax, k, r, n, il, ir;
  double x, y, z, tot;
  double *arr, *cnt;
  double gam, mag;

  hs = hist_open(hs0->rows, hs0->xmin,
      hs0->xmin + hs0->n * hs0->dx, hs0->dx);

  n = hs->n;
  xnew(cnt, n + 1);

  for ( r = 0; r < hs->rows; r++ ) {
    /* bootstrap for histogram r */
    arr = hs0->arr + r * hs->n;

    /* count the total for histogram r */
    tot = 0;
    imin = n;
    imax = 0;
    cnt[0] = 0;
    for ( i = 0; i < n; i++ ) {
      x = arr[i];
      cnt[i + 1] = cnt[i] + x;
      if ( x > 0 ) {
        if ( i < imin ) imin = i;
        if ( i + 1 > imax ) imax = i + 1;
      }
    }
    tot = cnt[n];

    /* bootstrapping based on Hub 2010 */
    x = randgaus();
    if ( tau != NULL && tau[r] > 0 ) {
      gam = exp(-1/tau[r]);
      mag = sqrt(1 - gam * gam);
    } else {
      gam = 0;
      mag = 1;
    }

    for ( i = 0, k = 0; k < tot; k++ ) {
      x = x * gam + mag * randgaus();
      /* convert the normally distributed variable to [0, 1] */
      y = 0.5 * (erf(x / sqrt(2)) + 1);

      /* translate y to the histogram frame */
      z = tot * y;
      /* find the bin i containing z using binary search
       * that is cnt[i] < z < cnt[i+1] */
      il = imin;
      ir = imax;
      while ( il < ir - 1 ) {
        i = (il + ir + 1) / 2;
        if ( z > cnt[i] ) {
          il = i;
        } else {
          ir = i;
        }
      }
      i = il;
      hs->arr[r * n + i] += 1;
    }
  }

  free(cnt);
  return hs;
}



/* reweight histograms by the autocorrelation time
 * `tau` is measured in terms of the number of frames */
static void reweightact(hist_t *hs, const double *tau)
{
  int i, j, n = hs->n;

  for ( i = 0; i < hs->rows; i++ ) {
    double w = tau[i] * 2 + 1;
    for ( j = 0; j < n; j++ ) {
      hs->arr[i*n + j] /= w;
    }
  }
}



/* compute the probability for selecting temperatures */
static int calcprob(hist_t *hs,
    const xdouble *beta, const xdouble *lnz,
    int id, double *proba)
{
  int i, n = hs->rows;
  double e, y, max;

  e = hs->xmin + (id + 0.5) * hs->dx;

  /* find the maximal of -beta * e - lnz */
  max = -DBL_MAX;
  for ( i = 0; i < n; i++ ) {
    y = (double) (-beta[i] * e - lnz[i]);
    if ( y > max ) max = y;
  }

  /* compute the cumulative distribution function */
  proba[0] = 0;
  for ( i = 0; i < n; i++ ) {
    y = (double) (-beta[i] * e - lnz[i] - max);
    proba[i+1] = proba[i] + exp(y);
  }

  return i;
}



/* select a temperature using the heat bath algorithm */
static int selectbeta(double *proba, int n)
{
  int i;
  double r;

  /* select the temperature */
  r = proba[n] * rand01();
  for ( i = 0; i < n; i++ ) {
    if ( r >= proba[i] && r < proba[i+1] ) {
      break;
    }
  }

  if ( i >= n ) {
    fprintf(stderr, "error! r %g, proba %g\n", r, proba[n]); getchar();
    return n - 1;
  }

  return i;
}



/* expanded ensemble bootstrapping with autocorrelation time `tau`
 * `tau` is measured in terms of the number of frames */
static hist_t *hist_eebootstrap(hist_t *hs0,
    const xdouble *beta, const xdouble *lnz,
    const double *tau)
{
  hist_t *hs;
  int K = hs0->rows;
  int i, j, r, n, il, ir;
  int *imin, *imax;
  double x;
  double *arr, **cnt;
  double *tot, *ctot;
  double **proba, *ntot;
  double *gam;

  hs = hist_open(K, hs0->xmin,
      hs0->xmin + hs0->n * hs0->dx, hs0->dx);

  n = hs->n;
  xnew(cnt, K);
  for ( r = 0; r < K; r++ ) {
    xnew(cnt[r], n + 1);
  }
  xnew(tot, K);
  xnew(ctot, K + 1);
  xnew(proba, n);
  for ( i = 0; i < n; i++ ) {
    xnew(proba[i], K + 1);
  }
  xnew(ntot, K);
  xnew(imin, K);
  xnew(imax, K);
  xnew(gam, K);

  /* compute tot, imin, imax, etc. */
  ctot[0] = 0;
  for ( r = 0; r < K; r++ ) {
    /* for histogram r */
    arr = hs0->arr + r * hs->n;

    /* count the total for histogram r */
    tot[r] = 0;
    imin[r] = n;
    imax[r] = 0;
    cnt[r][0] = 0;
    for ( i = 0; i < n; i++ ) {
      x = arr[i];
      cnt[r][i + 1] = cnt[r][i] + x;
      if ( x > 0 ) {
        if ( i < imin[r] ) imin[r] = i;
        if ( i + 1 > imax[r] ) imax[r] = i + 1;
      }
    }
    tot[r] = cnt[r][n];

    /* cumulative version of `tot` */
    ctot[r + 1] = ctot[r] + tot[r];

    gam[r] = (tau != NULL && tau[r] > 0) ? exp(-1/tau[r]) : 0;
  }

  for ( i = 0; i < n; i++ ) {
    calcprob(hs, beta, lnz, i, proba[i]);
  }

  /* bootstrapping */
  r = (int) (rand01() * K); /* temperature index */
  i = -1; /* frame index */
  for ( j = 0; j < ctot[K]; j++ ) {
    /* randomly pick a frame in the current temperature */
    if ( i < 0 || rand01() >= gam[r] ) {
      x = tot[r] * rand01();
      /* find the bin i containing x using binary search
       * that is cnt[r][i] < x < cnt[r][i+1] */
      il = imin[r];
      ir = imax[r];
      while ( il < ir - 1 ) {
        i = (il + ir + 1) / 2;
        if ( x > cnt[r][i] ) {
          il = i;
        } else {
          ir = i;
        }
      }
      i = il;
    }

    /* 2. try to change the temperature */
    r = selectbeta(proba[i], K);

    hs->arr[r * n + i] += 1;
    ntot[r] += 1;
  }

  for ( r = 0; r < hs->rows; r++ ) {
    fprintf(stderr, "%4d %10.7f %10.f %10.0f\n",
        r, (double) beta[r], tot[r], ntot[r]);
  }

  for ( r = 0; r < K; r++ ) {
    free(cnt[r]);
  }
  free(cnt);
  free(tot);
  free(ctot);
  for ( i = 0; i < n; i++ ) {
    free(proba[i]);
  }
  free(proba);
  free(ntot);
  free(imin);
  free(imax);
  return hs;
}



int main(int argc, char **argv)
{
  model_t m[1];
  hist_t *hs = NULL, *hs0 = NULL;
  int i, nbeta;
  unsigned flags = 0;
  xdouble *beta, *lnz;
  double *tcorr = NULL; /* array of correlation time */


  model_default(m);
  model_doargs(m, argc, argv);
  if ( m->rmcom ) flags |= WHAM_RMCOM;

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  if ( m->loadprev ) {
    /* load from the existing histogram */
    if ( (hs = hist_initf(m->fnhis, HIST_INT | HIST_VERBOSE)) == NULL ) {
      return -1;
    }

    getls(m->fninp, &nbeta, &beta);
    if ( nbeta != hs->rows ) {
      fprintf(stderr, "%s: different rows %d vs %d\n",
          m->fnhis, nbeta, hs->rows);
      hist_close(hs);
      return -1;
    }
  } else {
    hs = mkhist(m->fninp, &beta, m->de, m->fnhis, m->radd,
                m->actmax, m->acmin, m->fnact);
    if ( hs == NULL ) {
      return -1;
    }
  }

  xnew(tcorr, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) tcorr[i] = 0;

  if ( m->bootstrap || m->eebootstrap || m->weightact ) {
    /* load autocorrelation time */
    if ( loadact(m->fnact, hs->rows, tcorr) != 0 ) {
      if ( m->fnact == NULL ) {
        fprintf(stderr, "no file for autocorrelation time\n");
      } else {
        fprintf(stderr, "cannot load autocorrelation time from %s\n",
            m->fnact);
      }
      free(tcorr);
      tcorr = NULL;
    }
  }

  if ( m->bootstrap ) {
    hs0 = hs;

    /* scramble the random number seed */
    mtscramble( time(NULL) + clock() );
    /* bootstrapping */
    hs = hist_bootstrap(hs0, tcorr);
    hist_close(hs0);
  }

  xnew(lnz, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) {
    lnz[i] = 0;
  }

  if ( m->weightact ) {
    reweightact(hs, tcorr);
  }

  whamx(hs, beta, lnz, flags, NULL,
      m->damp, m->mdiis_nbases,
      m->mdiis_update_method, m->mdiis_threshold,
      m->itmin, m->itmax, m->tol, m->verbose,
      m->fnlndos, m->fneav, m->wham_method);

  if ( m->eebootstrap ) {
    hs0 = hs;

    /* scramble the random number seed */
    mtscramble( time(NULL) + clock() );
    /* expanded ensemble bootstrapping
     * if the histogram has already been reweighted by `weightact`
     * do not apply the correlation time */
    hs = hist_eebootstrap(hs0, beta, lnz,
        m->weightact ? NULL : tcorr);
    hist_close(hs0);

    /* do WHAM for the bootstrap sample */
    whamx(hs, beta, lnz, flags, NULL,
        m->damp, m->mdiis_nbases,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmin, m->itmax, m->tol, m->verbose,
        m->fnlndos, m->fneav, m->wham_method);
  }

  if ( m->verbose ) {
    for ( i = 0; i < hs->rows; i++ ) {
      double tot, eav, var;
      eav = hist_getave(hs, i, &tot, &var);
      printf("%3d %10.7f %14.7f %8.0f %15.7f %14.7f\n",
          i, (double) beta[i], (double) lnz[i], tot, eav, sqrt(var));
    }
  }

  hist_close(hs);
  free(beta);
  free(lnz);
  free(tcorr);
  return 0;
}

