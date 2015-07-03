/* WHAM for a set of xvg files */
#define ENABLE_MDIIS
#include "../wham2_xdbl.h"
#include <time.h>
#include "../mtrand.h"
#define MTRAND
#include "xvg.h"
#define WHAM
#include "../whammodel.h"
#include "ls2util.h"




#define EVCNT  2



/* construct the histogram */
static hist2_t *mkhist2(const char *fnls,
    xdouble **beta, xdouble **bpres, double de, double dv,
    const char *fnhis, double radd)
{
  hist2_t *hs;
  int i, j, nbp, k;
  char **fns;
  xvg_t **xvg = NULL;
  double evmin[EVCNT] = {1e30, 1e30}, evmax[EVCNT] = {-1e30, -1e30};
  double evmin1[EVCNT], evmax1[EVCNT];

  /* scramble the random number seed */
  mtscramble( time(NULL) );

  if ( (fns = getls(fnls, &nbp, beta, bpres)) == NULL ) {
    return NULL;
  }

  xnew(xvg, nbp);
  for ( i = 0; i < nbp; i++ ) {
    if ( (xvg[i] = xvg_load(fns[i], radd)) == NULL ) {
      exit(1);
    }
    xvg_minmax(xvg[i], evmin1, evmax1);
    for ( k = 0; k < EVCNT; k++ ) {
      if ( evmin1[k] < evmin[k] ) {
        evmin[k] = evmin1[k];
      }
      if ( evmax1[k] > evmax[k] ) {
        evmax[k] = evmax1[k];
      }
    }
  }

  evmin[0] = ((int) (evmin[0] / de) - 1) * de;
  evmax[0] = ((int) (evmax[0] / de) + 1) * de;
  evmin[1] = ((int) (evmin[1] / dv) - 1) * dv;
  evmax[1] = ((int) (evmax[1] / dv) + 1) * dv;

  hs = hist2_open(nbp, evmin[0], evmax[0], de,
      evmin[1], evmax[1], dv);

  for ( i = 0; i < nbp; i++ ) {
    for ( j = 0; j < xvg[i]->n; j++ ) {
      hist2_add1(hs, i,
          xvg[i]->y[0][j], xvg[i]->y[2][j], 1.0, HIST2_VERBOSE);
    }
  }

  hist2_save(hs, fnhis, HIST2_ADDAHALF | HIST2_NOZEROES | HIST2_VERBOSE);

  for ( i = 0; i < nbp; i++ ) {
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
  hist2_t *hs = NULL;
  int i, nbp;
  xdouble *beta, *bpres, *lnz;

  model_default(m);
  m->de = 10.0;
  model_doargs(m, argc, argv);

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  /* try to load from the existing histogram */
  if ( m->loadprev ) {
    hs = hist2_initf(m->fnhis2, HIST2_INT | HIST2_VERBOSE);
    if ( hs == NULL ) {
      return -1;
    }

    getls(m->fninp, &nbp, &beta, &bpres);
    if ( nbp != hs->rows ) {
      fprintf(stderr, "%s: different rows %d vs %d\n",
          m->fnhis2, nbp, hs->rows);
      hist2_close(hs);
      return -1;
    }
  } else {
    hs = mkhist2(m->fninp, &beta, &bpres, m->de, m->dv,
        m->fnhis2, m->radd);
    if ( hs == NULL ) {
      return -1;
    }
  }

  xnew(lnz, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) {
    lnz[i] = 0;
  }

  wham2x(hs, beta, bpres, lnz,
      m->damp, m->mdiis_nbases,
      m->mdiis_update_method, m->mdiis_threshold,
      m->itmin, m->itmax, m->tol, m->verbose,
      m->fnlndos2, m->fneav2, m->wham_method);

  if ( m->verbose ) {
    for ( i = 0; i < hs->rows; i++ ) {
      double tot, eav, vav, see, sev, svv, sige, sigv;
      eav = hist2_getave(hs, i, &tot, &vav, &see, &sev, &svv);
      sige = sqrt(see);
      sigv = sqrt(svv);
      printf("%3d %10.7f %10.7f %14.7f %8.0f "
          "%15.7f %14.7f %15.7f %14.7f %15.7f\n",
          i, (double) beta[i], (double) (bpres[i]/beta[i]),
          (double) lnz[i], tot,
          eav, sige, vav, sigv, sev);
    }
  }

  hist2_close(hs);
  free(beta);
  free(lnz);
  return 0;
}

