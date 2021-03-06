/* WHAM for two-dimensional Ising model */
#define ENABLE_MDIIS
#include "../wham.h"
#ifndef IS2_LB
#define IS2_LB 6
#endif
#include "is2.h"
#include <time.h>
#define WHAM
#define IS2_MODEL
#include "../whammodel.h"



double xmin = -2*IS2_N - 2, xmax = 2;



static hist_t *is2_simul(model_t *m, double *beta)
{
  int iT, id, h, istep;
  is2_t **is;
  hist_t *hs;

  hs = hist_open(m->nT, xmin, xmax, m->de);

  xnew(is, m->nT);
  for ( iT = 0; iT < m->nT; iT++ ) {
    is[iT] = is2_open(IS2_L);
    IS2_SETPROBA(is[iT], beta[iT]);
  }

  /* randomize the initial state */
  mtscramble( time(NULL) );

  /* do the simulations */
  for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
    /* MC for each temperature */
    for ( iT = 0; iT < m->nT; iT++ ) {
      IS2_PICK(is[iT], id, h);
      if ( h < 0 || mtrand() <= is[iT]->uproba[h] ) {
        IS2_FLIP(is[iT], id, h);
      }
    }

    if ( m->re ) {
      /* replica exchange: randomly swap configurations of
       * two neighboring temperatures */
      int jT, acc;
      double dbdE, r;
      is2_t *istmp;
      unsigned utmp;

      iT = (int) (rand01() * (m->nT - 1));
      jT = iT + 1;
      dbdE = (beta[iT] - beta[jT]) * (is[iT]->E - is[jT]->E);
      acc = 0;
      if ( dbdE >= 0 ) {
        acc = 1;
      } else {
        r = rand01();
        if ( r < exp(dbdE) ) {
          acc = 1;
        }
      }
      if ( acc ) {
        /* swap the models */
        istmp = is[iT], is[iT] = is[jT], is[jT] = istmp;
        /* swap the transition probabilities */
        utmp = is[iT]->uproba[2], is[iT]->uproba[2] = is[jT]->uproba[2], is[jT]->uproba[2] = utmp;
        utmp = is[iT]->uproba[4], is[iT]->uproba[4] = is[jT]->uproba[4], is[jT]->uproba[4] = utmp;
      }
    }

    if ( istep <= m->nequil ) continue;
    //for ( iT = 0; iT < m->nT; iT++ ) printf("iT %d, ep %d\n", iT, is[iT]->E);
    //printf("hs->xmin %g\n", hs->xmin); getchar();
    for ( iT = 0; iT < m->nT; iT++ ) {
      hist_add1(hs, iT, is[iT]->E, 1.0, HIST_VERBOSE);
    }
  }

  hist_save(hs, m->fnhis, HIST_ADDAHALF | HIST_NOZEROES | HIST_VERBOSE);
  fprintf(stderr, "simulation ended in %d steps, doing WHAM\n", m->nsteps);

  for ( iT = 0; iT < m->nT; iT++ ) {
    is2_close( is[iT] );
  }
  return hs;
}



static void model_default_is2(model_t *m)
{
  model_default(m);
  m->de = 4;
  m->nT = 80;
  m->Tmin = 1.5;
  m->Tdel = 0.02;
  m->nequil = 100000;
  m->nsteps = 10000000;
}



int main(int argc, char **argv)
{
  model_t m[1];
  hist_t *hs;
  double *beta, *lnz, *lnzref, *lnz1;
  int iT;
  unsigned flags = 0;

  model_default_is2(m);
  model_doargs(m, argc, argv);
  if ( m->rmcom ) flags |= WHAM_RMCOM;

  xnew(beta, m->nT);
  xnew(lnz, m->nT);
  xnew(lnzref, m->nT);
  xnew(lnz1, m->nT);
  for ( iT = 0; iT < m->nT; iT++ ) {
    beta[iT] = 1./(m->Tmin + m->Tdel * iT);
    lnz[iT] = 0;
    lnz1[iT] = 0;
    lnzref[iT] = is2_exact(IS2_L, IS2_L, beta[iT], NULL, NULL);
  }

  if ( m->loadprev ) {
    /* load from existing histogram */
    if ( (hs = hist_initf(m->fnhis, HIST_INT | HIST_VERBOSE)) == NULL ) {
      return -1;
    }
  } else {
    hs = is2_simul(m, beta);
  }

  whamx(hs, beta, lnz, flags, NULL,
      m->damp, m->mdiis_nbases,
      m->mdiis_update_method, m->mdiis_threshold,
      m->itmin, m->itmax, m->tol, m->verbose,
      m->fnlndos, m->fneav, m->wham_method);

  if ( m->verbose ) {
    for ( iT = 0; iT < m->nT; iT++ ) {
      printf("%3d %10.7f %14.7f %14.7f\n",
         iT, beta[iT], lnz[iT], lnzref[iT] - lnzref[0]);
    }
  }

  if ( m->debug ) {
    /* increase the verbose level, using the converged `lnz`
     * as the reference to compute the error */
    whamx(hs, beta, lnz1, flags, lnz,
        m->damp, m->mdiis_nbases,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmin, m->itmax, m->tol, m->verbose + 1,
        m->fnlndos, m->fneav, m->wham_method);
  }

  hist_close(hs);
  free(beta);
  free(lnz);
  free(lnz1);
  free(lnzref);
  return 0;
}

