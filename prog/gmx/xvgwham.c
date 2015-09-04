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
  int i;
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

  fgets(s, sizeof s, fp);
  sscanf(s + 1, "%d%lf", &i, &dx);

  for ( i = 0; i < nbet; i++ ) {
    fgets(s, sizeof s, fp);
    sscanf(s, "%lf", &tcorr[i]);
    tcorr[i] /= dx;
    //fprintf(stderr, "%d: %g %g\n", i, tcorr[i], tcorr[i] * dx);
  }

  fclose(fp);
  return 0;
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



/* bootstrapping with autocorrelation time `tau`
 * `tau` is measured in terms of the number of frames */
static hist_t *hist_bootstrap(hist_t *hs0, double *tau)
{
  hist_t *hs;
  int i, imin, imax, k, r, n, il, ir;
  double x, tot, *arr, *cnt;
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

    /* bootstrapping */
    gam = (tau != NULL && tau[r] > 0 ? exp(-1/tau[r]) : 0);
    for ( i = 0, k = 0; k < tot; k++ ) {
      /* selectively update the current frame */
      if ( rand01() > gam || k == 0 ) {
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



int main(int argc, char **argv)
{
  model_t m[1];
  hist_t *hs = NULL;
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

  if ( m->bootstrap ) {
    hist_t *hs0 = hs;

    /* load autocorrelation time */
    xnew(tcorr, hs->rows);
    if ( loadact(m->fnact, hs->rows, tcorr) != 0 ) {
      fprintf(stderr, "cannot load autocorrelation time from %s\n", m->fnact);
      free(tcorr);
      tcorr = NULL;
    }

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

  whamx(hs, beta, lnz, flags, NULL,
      m->damp, m->mdiis_nbases,
      m->mdiis_update_method, m->mdiis_threshold,
      m->itmin, m->itmax, m->tol, m->verbose,
      m->fnlndos, m->fneav, m->wham_method);

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

