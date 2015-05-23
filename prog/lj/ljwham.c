/* test program for WHAM */
#define LJ_MODEL
#include "../whammodel.h"
#define WHAM_MDIIS
#include "../wham.h"
#include "lj.h"
#include <time.h>




static void model_default_lj(model_t *m)
{
  model_default(m);
  m->np = 108;
  m->rho = 0.3;
  m->rcdef = 2.5;
  m->mddt = 0.002;
  m->thdt = 0.02;
  m->pdt = 1e-5;
  m->de = 0.1;
  m->emin = -6.0;
  m->emax = 1.0;
  m->nT = 10;
  m->Tmin = 0.7;
  m->Tdel = 0.1;
  m->nequil = 5000;
  m->nsteps = 100000;
}



int main(int argc, char **argv)
{
  model_t m[1];
  hist_t *hs;
  double *beta, *lnz;
  int iT;

  model_default_lj(m);
  model_doargs(m, argc, argv);

  xnew(beta, m->nT);
  xnew(lnz, m->nT);
  for ( iT = 0; iT < m->nT; iT++ ) {
    beta[iT] = 1./(m->Tmin + m->Tdel * iT);
    lnz[iT] = 0;
  }

  if ( m->loadprev ) {
    /* load from existing histogram */
    if ( (hs = hist_initf(m->fnhis)) == NULL ) {
      return -1;
    }
  } else {
    int istep;
    lj_t **lj;

    hs = hist_open(m->nT, m->np * m->emin, m->np * m->emax, m->de);

    /* randomize the initial state */
    mtscramble( time(NULL) );

    xnew(lj, m->nT);
    for ( iT = 0; iT < m->nT; iT++ ) {
      lj[iT] = lj_open(m->np, m->rho, m->rcdef);
    }

    /* do simulations */
    for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
      for ( iT = 0; iT < m->nT; iT++ ) {
        lj_vv(lj[iT], m->mddt);
        lj[iT]->ekin = lj_vrescale(lj[iT], 1/beta[iT], m->thdt);
        hist_add1(hs, iT, lj[iT]->epot, 1.0, HIST_VERBOSE);
      }
      if ( istep <= m->nequil ) continue;
    }

    hist_save(hs, m->fnhis, HIST_ADDAHALF);
    fprintf(stderr, "simulation ended, doing WHAM\n");

    for ( iT = 0; iT < m->nT; iT++ ) {
      lj_close( lj[iT] );
    }
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    wham(hs, beta, lnz,
        m->itmax, m->tol, m->verbose, m->fnlndos, m->fneav);
  } else {
    wham_mdiis(hs, beta, lnz, m->mdiis_nbases, m->mdiis_damp,
        m->itmax, m->tol, m->verbose, m->fnlndos, m->fneav);
  }

  hist_close(hs);
  free(beta);
  free(lnz);
  return 0;
}

