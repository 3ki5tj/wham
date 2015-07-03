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





/* construct the histogram */
static hist_t *mkhist(const char *fnls,
    xdouble **beta, double de,
    const char *fnhis, double radd)
{
  hist_t *hs;
  int i, j, nbeta;
  char **fns;
  xvg_t **xvg = NULL;
  double xmin[2], xmax[2], emin = 1e30, emax = -1e30;

  /* scramble the random number seed */
  mtscramble( time(NULL) );

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

  for ( i = 0; i < nbeta; i++ ) {
    xvg_close(xvg[i]);
    free(fns[i]);
  }
  free(xvg);
  free(fns);

  return hs;
}



/* bootstrapping */
static hist_t *hist_bootstrap(hist_t *hs0)
{
  hist_t *hs;
  int i, imin, imax, k, r, n, il, ir;
  double x, tot, *arr, *cnt;

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
    for ( k = 0; k < tot; k++ ) {
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
  xdouble *beta, *lnz;

  model_default(m);
  model_doargs(m, argc, argv);

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
    hs = mkhist(m->fninp, &beta, m->de, m->fnhis, m->radd);
    if ( hs == NULL ) {
      return -1;
    }
  }

  if ( m->bootstrap ) {
    hist_t *hs0 = hs;

    /* scramble the random number seed */
    mtscramble( time(NULL) );
    /* bootstrapping */
    hs = hist_bootstrap(hs0);
    hist_close(hs0);
  }

  xnew(lnz, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) {
    lnz[i] = 0;
  }

  whamx(hs, beta, lnz,
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
  return 0;
}

