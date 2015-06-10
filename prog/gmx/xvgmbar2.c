/* MBAR for a set of xvg files */
#include <time.h>
#include "../mtrand.h"
#define MTRAND
#include "xvg.h"
#define ENABLE_MDIIS
#include "mbar2.h"
#define MBAR
#include "../whammodel.h"
#include "ls2util.h"





/* construct the histogram */
static xvg_t **mkxvg(const char *fnls,
    int *nbp, double **beta, double **bpres,
    double radd)
{
  int i;
  char **fns;
  xvg_t **xvg = NULL;

  /* scramble the random number seed */
  mtscramble( time(NULL) );

  if ( (fns = getls(fnls, nbp, beta, bpres)) == NULL ) {
    return NULL;
  }

  xnew(xvg, *nbp);
  for ( i = 0; i < *nbp; i++ ) {
    if ( (xvg[i] = xvg_load(fns[i], radd)) == NULL ) {
      exit(1);
    }
  }

  for ( i = 0; i < *nbp; i++ ) {
    free(fns[i]);
  }
  free(fns);

  return xvg;
}



int main(int argc, char **argv)
{
  model_t m[1];
  xvg_t **xvg = NULL;
  int i, nbp;
  double *beta, *bpres, *lnz;

  model_default(m);
  model_doargs(m, argc, argv);

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  xvg = mkxvg(m->fninp, &nbp, &beta, &bpres, m->radd);
  if ( xvg == NULL ) {
    return -1;
  }
  fprintf(stderr, "trajectories loaded, doing MBAR...\n");

  xnew(lnz, nbp);
  for ( i = 0; i < nbp; i++ ) {
    lnz[i] = 0;
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    mbar2(nbp, xvg, beta, bpres, lnz,
        m->itmax, m->tol, m->itmin, m->verbose);
  } else {
    mbar2_mdiis(nbp, xvg, beta, bpres, lnz,
        m->mdiis_nbases, m->mdiis_damp,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmax, m->tol, m->itmin, m->verbose);
  }

  if ( m->verbose ) {
    for ( i = 0; i < nbp; i++ ) {
      printf("%3d %10.7f %10.7f %14.7f %8d\n",
          i, beta[i], bpres[i]/beta[i], lnz[i], xvg[i]->n);
    }
  }

  for ( i = 0; i < nbp; i++ ) {
    free(xvg[i]);
  }
  free(xvg);
  free(beta);
  free(lnz);
  return 0;
}

