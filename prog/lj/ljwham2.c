/* test program for the two-dimensional WHAM */
#define LJ_MODEL
#include "../whammodel.h"
#define WHAM2_MDIIS
#include "../wham2.h"
#include "lj.h"
#include <time.h>




static void model_default_lj2(model_t *m)
{
  model_default(m);
  m->np = 108;
  m->rho = 0.7;
  m->rcdef = 2.5;
  m->mddt = 0.002;
  m->thdt = 0.02;
  m->pdt = 1e-5;
  m->de = 0.5;
  m->emin = -8.0;
  m->emax = 2.0;
  m->dv = 0.2;
  m->vmin = 0.0;
  m->vmax = 4.0;
  m->nT = 5;
  m->Tmin = 1.0;
  m->Tdel = 0.1;
  m->nP = 5;
  m->Pmin = 1.0;
  m->Pdel = 0.2;
  m->nequil = 5000;
  m->nsteps = 100000;
}



int main(int argc, char **argv)
{
  model_t m[1];
  hist2_t *hs;
  double *beta, *pres, *bp, *lnz;
  int it, ip, itp, ntp;

  model_default_lj2(m);
  model_doargs(m, argc, argv);

  ntp = m->nT * m->nP;
  xnew(beta, ntp);
  xnew(pres, ntp);
  xnew(bp, ntp);
  xnew(lnz, ntp);

  /* initialize the systems */
  itp = 0;
  for ( it = 0; it < m->nT; it++ ) {
    double T = m->Tmin + m->Tdel * it;
    for ( ip = 0; ip < m->nP; ip++ ) {
      double P = m->Pmin + m->Pdel * ip;
      beta[itp] = 1./T;
      pres[itp] = P;
      bp[itp] = beta[itp] * pres[itp];
      lnz[itp] = 0;
      itp++;
    }
  }

  if ( m->loadprev ) {
    /* load from existing histogram */
    if ( (hs = hist2_initf(m->fnhis2)) == NULL ) {
      return -1;
    }
  } else {
    int istep;
    lj_t **lj;

    hs = hist2_open(ntp, m->np * m->emin, m->np * m->emax, m->de,
        m->np * m->vmin, m->np * m->vmax, m->dv);

    /* randomize the initial state */
    mtscramble( time(NULL) );

    xnew(lj, ntp);
    for ( itp = 0; itp < ntp; itp++ ) {
      lj[itp] = lj_open(m->np, m->rho, m->rcdef);
    }

    /* do simulations */
    for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
      for ( itp = 0; itp < ntp; itp++ ) { /* loop over replicas */
        lj_vv(lj[itp], m->mddt);
        lj[itp]->ekin = lj_vrescale(lj[itp], 1/beta[itp], m->thdt);
        if ( istep % 5 == 0 ) {
          lj_langp0(lj[itp], m->pdt, 1/beta[itp], pres[itp], 0);
        }
      }
      if ( istep <= m->nequil ) continue;

      for ( itp = 0; itp < ntp; itp++ ) {
        hist2_add1(hs, itp, lj[itp]->epot, lj[itp]->vol,
            1.0, HIST2_VERBOSE);
      }
    }

    hist2_save(hs, m->fnhis2, HIST2_ADDAHALF);
    fprintf(stderr, "simulation ended, doing WHAM\n");

    for ( itp = 0; itp < ntp; itp++ ) {
      lj_close( lj[itp] );
    }
  }

  /* do WHAM */
  if ( m->wham_method == WHAM_DIRECT ) {
    wham2(hs, beta, bp, lnz,
        m->itmax, m->tol, m->verbose, m->fnlndos2, m->fneav2);
  } else {
    wham2_mdiis(hs, beta, bp, lnz, m->mdiis_nbases, m->mdiis_damp,
        m->itmax, m->tol, m->verbose, m->fnlndos2, m->fneav2);
  }

  /* clean up */
  hist2_close(hs);
  free(beta);
  free(pres);
  free(bp);
  free(lnz);
  return 0;
}

