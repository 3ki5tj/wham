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





/* load trajectories */
static xvg_t **mkxvg(const char *fnls,
    int *nbeta, double **beta, double radd)
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



/* estimate free energies */
static int estimate(int nbeta, xvg_t **xvg,
    const double *beta)
{
  double *eav, *var, *skw;
  double *lnza, *lnzb, *lnzc;
  double dlnza, dlnzb, dlnzc, e, db;
  int i, j, n;

  xnew(eav, nbeta);
  xnew(var, nbeta);
  xnew(skw, nbeta);
  /* compute the variances */
  for ( j = 0; j < nbeta; j++ ) {
    eav[j] = var[j] = skw[j] = 0;
    n = xvg[j]->n;

    /* compute the mean */
    for ( i = 0; i < n; i++ ) {
      e = xvg[j]->y[0][i];
      eav[j] += e;
    }
    eav[j] /= n;

    /* compute the variance and skew */
    for ( i = 0; i < n; i++ ) {
      e = xvg[j]->y[0][i] - eav[j];
      var[j] += e * e;
      skw[j] += e * e * e;
    }
    if ( n > 1 ) {
      var[j] /= n - 1;
    }
    if ( n > 2 ) {
      skw[j] *= 1.0 * n / (n - 1) / (n - 2);
    }
  }

  xnew(lnza, nbeta);
  xnew(lnzb, nbeta);
  xnew(lnzc, nbeta);
  lnza[0] = lnzb[0] = lnzc[0] = 0;
  for ( j = 0; j < nbeta - 1; j++ ) {
    db = beta[j+1] - beta[j];
    dlnza = -(eav[j+1] + eav[j]) / 2 * db;
    dlnzb = dlnza - (var[j+1] - var[j]) / 12 * db * db;
    dlnzc = dlnza + (skw[j+1] + skw[j]) / 24 * db * db * db;
    lnza[j+1] = lnza[j] + dlnza;
    lnzb[j+1] = lnzb[j] + dlnzb;
    lnzc[j+1] = lnzc[j] + dlnzc;
  }

  for ( j = 0; j < nbeta; j++ ) {
    printf("%3d %10.7f %14.7f %14.7f %14.7f %15.7f %14.7f %8d\n",
        j, beta[j], lnza[j], lnzb[j], lnzc[j],
        eav[j], sqrt(var[j]), xvg[j]->n);
  }

  free(eav);
  free(var);
  free(skw);
  free(lnza);
  free(lnzb);
  free(lnzc);
  return 0;
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

  xvg = mkxvg(m->fninp, &nbeta, &beta, m->radd);
  if ( xvg == NULL ) {
    return -1;
  }

  if ( m->bootstrap ) {
    xvg_t *xvg0;

    /* scramble the random number seed */
    mtscramble( time(NULL) );
    /* bootstrapping */
    for ( i = 0; i < nbeta; i++ ) {
      xvg0 = xvg[i];
      xvg[i] = xvg_bootstrap( xvg0 );
      xvg_close(xvg0);
    }
  }

  fprintf(stderr, "trajectories loaded, doing MBAR...\n");

  xnew(lnz, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    lnz[i] = 0;
  }

  if ( m->estimate ) {
    /* estimate */
    estimate(nbeta, xvg, beta);
  } else {
    /* do MBAR */
    mbarx(nbeta, xvg, beta, lnz,
        m->damp, m->mdiis_nbases,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmin, m->itmax, m->tol, m->verbose, m->mbar_method);

    if ( m->verbose ) {
      for ( i = 0; i < nbeta; i++ ) {
        printf("%3d %10.7f %14.7f %8d\n",
            i, beta[i], lnz[i], xvg[i]->n);
      }
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

