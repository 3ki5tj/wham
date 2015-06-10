/* test program for WHAM */
#define ENABLE_MDIIS
#include "../wham.h"
#include "lj.h"
#include <time.h>
#define LJ_MODEL
#include "../whammodel.h"




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
  m->nsteps = 500000;
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
    if ( (hs = hist_initf(m->fnhis, HIST_INT | HIST_VERBOSE)) == NULL ) {
      return -1;
    }
  } else {
    double emin, emax;
    int istep;
    lj_t **lj;

    /* round emin and emax to multiples of m->de */
    emin = (int) (m->np * m->emin / m->de) * m->de;
    emax = (int) (m->np * m->emax / m->de) * m->de;

    hs = hist_open(m->nT, emin, emax, m->de);

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
      }

      if ( m->re ) {
        /* replica exchange: randomly swap configurations of
         * two neighboring temperatures */
        int jT, acc;
        double dbdE, r;
        lj_t *ljtmp;

        iT = (int) (rand01() * (m->nT - 1));
        jT = iT + 1;
        dbdE = (beta[iT] - beta[jT]) * (lj[iT]->epot - lj[jT]->epot);
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
          double scl = sqrt( beta[iT]/beta[jT] ); /* sqrt(Tj/Ti) */
          int i;

          /* scale the velocities */
          for ( i = 0; i < m->np; i++ ) {
            vsmul(lj[iT]->v[i], scl);
            vsmul(lj[jT]->v[i], 1/scl);
          }
          /* swap the models */
          ljtmp = lj[iT], lj[iT] = lj[jT], lj[jT] = ljtmp;
        }
      }

      if ( istep <= m->nequil ) continue;

      for ( iT = 0; iT < m->nT; iT++ ) {
        hist_add1(hs, iT, lj[iT]->epot, 1.0, HIST_VERBOSE);
      }
    }

    hist_save(hs, m->fnhis, HIST_ADDAHALF);
    fprintf(stderr, "simulation ended, doing WHAM\n");

    for ( iT = 0; iT < m->nT; iT++ ) {
      lj_close( lj[iT] );
    }
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    wham(hs, beta, lnz,
        m->itmax, m->tol, m->itmin, m->verbose,
        m->fnlndos, m->fneav);
  } else {
    wham_mdiis(hs, beta, lnz,
        m->mdiis_nbases, m->mdiis_damp,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmax, m->tol, m->itmin, m->verbose,
        m->fnlndos, m->fneav);
  }

  hist_close(hs);
  free(beta);
  free(lnz);
  return 0;
}

