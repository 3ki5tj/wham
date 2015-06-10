/* MBAR for a set of xvg files */
#include <time.h>
#include "../mtrand.h"
#define MTRAND
#include "xvg.h"
#define ENABLE_MDIIS
#include "mbar.h"
#define MBAR
#include "../whammodel.h"
#include "lsutil.h"





/* construct the histogram */
static xvg_t **mkxvg(const char *fnls,
    int *nbeta, double **beta, double de,
    double radd)
{
  int i;
  char **fns;
  xvg_t **xvg = NULL;

  /* scramble the random number seed */
  mtscramble( time(NULL) );

  if ( (fns = getls(fnls, nbeta, beta)) == NULL ) {
    return NULL;
  }

  xnew(xvg, *nbeta);
  for ( i = 0; i < *nbeta; i++ ) {
    if ( (xvg[i] = xvg_load(fns[i], radd)) == NULL ) {
      exit(1);
    }
  }

  for ( i = 0; i < *nbeta; i++ ) {
    free(fns[i]);
  }
  free(fns);

  return xvg;
}



int main(int argc, char **argv)
{
  model_t m[1];
  xvg_t **xvg = NULL;
  int i, nbeta;
  double *beta, *lnz;

  model_default(m);
  model_doargs(m, argc, argv);

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  xvg = mkxvg(m->fninp, &nbeta, &beta, m->de, m->radd);
  if ( xvg == NULL ) {
    return -1;
  }
  fprintf(stderr, "trajectories loaded, doing MBAR...\n");

  xnew(lnz, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    lnz[i] = 0;
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    mbar(nbeta, xvg, beta, lnz,
        m->itmax, m->tol, m->itmin, m->verbose);
  } else {
    mbar_mdiis(nbeta, xvg, beta, lnz,
        m->mdiis_nbases, m->mdiis_damp,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmax, m->tol, m->itmin, m->verbose);
  }

  if ( m->verbose ) {
    for ( i = 0; i < nbeta; i++ ) {
      printf("%3d %10.7f %14.7f %8d\n",
          i, beta[i], lnz[i], xvg[i]->n);
    }
  }

  for ( i = 0; i < nbeta; i++ ) {
    free(xvg[i]);
  }
  free(xvg);
  free(beta);
  free(lnz);
  return 0;
}

