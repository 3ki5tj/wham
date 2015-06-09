/* WHAM for a set of xvg files */
#define ENABLE_MDIIS
#include "../wham.h"
#include <time.h>
#include "../mtrand.h"
#define MTRAND
#include "xvg.h"
#include "../whammodel.h"
#include "lsutil.h"





/* construct the histogram */
static hist_t *mkhist(const char *fnls,
    double **beta, double de,
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



int main(int argc, char **argv)
{
  model_t m[1];
  hist_t *hs = NULL;
  int i, nbeta;
  double *beta, *lnz;

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

  xnew(lnz, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) {
    lnz[i] = 0;
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    wham(hs, beta, lnz,
        m->itmax, m->tol, m->verbose, m->fnlndos, m->fneav);
  } else {
    wham_mdiis(hs, beta, lnz,
        m->mdiis_nbases, m->mdiis_damp,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmax, m->tol, m->verbose, m->fnlndos, m->fneav);
  }

  if ( m->verbose ) {
    for ( i = 0; i < hs->rows; i++ ) {
      double tot, eav, var;
      eav = hist_getave(hs, i, &tot, &var);
      fprintf(stderr, "%3d %8.5f %10.3f %8.0f %11.4f(%10.4f)\n",
          i, beta[i], lnz[i], tot, eav, sqrt(var));
    }
  }

  hist_close(hs);
  free(beta);
  free(lnz);
  return 0;
}

