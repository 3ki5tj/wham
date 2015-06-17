/* test program for WHAM */
#define ENABLE_MDIIS
#include "../wham.h"
#include "lj.h"
#include "ljeos.h"
#include <time.h>
#define WHAM
#define LJ_MODEL
#include "../whammodel.h"




static void model_default_lj(model_t *m)
{
  model_default(m);
  m->nn = 108;
  m->rho = 0.7;
  m->rcdef = 2.5;
  m->mcamp = 0.5;
  m->mddt = 0.002;
  m->thdt = 0.02;
  m->pdt = 1e-5;
  m->de = 0.1;
  m->emin = -6.0;
  m->emax = 1.0;
  m->nT = 10;
  m->Tmin = 0.7;
  m->Tdel = 0.1;
  m->nstadj = 5000;
  m->nequil = 25000;
  m->nsteps = 1000000;
  m->simul = SIMUL_MC;
}



static hist_t *lj_domc(model_t *m, const double *beta)
{
  double emin, emax;
  int istep, iT, jT, acc;
  double *mctot, *mcacc, *mcamp;
  double retot = DBL_MIN, reacc = 0;
  lj_t **lj, *ljtmp;
  hist_t *hs;

  /* round emin and emax to multiples of m->de */
  emin = (int) (m->nn * m->emin / m->de) * m->de;
  emax = (int) (m->nn * m->emax / m->de) * m->de;

  hs = hist_open(m->nT, emin, emax, m->de);

  /* randomize the initial state */
  mtscramble( time(NULL) );

  xnew(lj, m->nT);
  for ( iT = 0; iT < m->nT; iT++ ) {
    lj[iT] = lj_open(m->nn, m->rho, m->rcdef);
    lj[iT]->dof = m->nn * D;
    lj_energy( lj[iT] );
  }

  xnew(mctot, m->nT);
  xnew(mcacc, m->nT);
  xnew(mcamp, m->nT);
  for ( iT = 0; iT < m->nT; iT++ ) {
    mctot[iT] = DBL_MIN;
    mcacc[iT] = 0;
    mcamp[iT] = m->mcamp;
  }

  /* do simulations */
  for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
    for ( iT = 0; iT < m->nT; iT++ ) {
      acc = lj_metro(lj[iT], mcamp[iT], beta[iT]);
      mcacc[iT] += acc;
      mctot[iT] += 1;
    }

    if ( m->re ) {
      /* replica exchange: randomly swap configurations of
       * two neighboring temperatures */
      double dbdE, r;

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
        /* swap the models */
        ljtmp = lj[iT], lj[iT] = lj[jT], lj[jT] = ljtmp;
      }

      if ( istep > m->nequil ) {
        reacc += acc;
        retot += 1;
      }
    }

    if ( istep <= m->nequil ) {
      /* adjust the MC move amplitude */
      if ( istep % m->nstadj == 0 ) {
        for ( iT = 0; iT < m->nT; iT++ ) {
          fprintf(stderr, "T %g, mcamp %g, mcacc %g%%\n",
              1/beta[iT], mcamp[iT], 100*mcacc[iT]/mctot[iT]);
          mcamp[iT] *= sqrt( mcacc[iT] / mctot[iT] / 0.5 );
        }
      }
      continue;
    }

    for ( iT = 0; iT < m->nT; iT++ ) {
      hist_add1(hs, iT, lj[iT]->epot, 1.0, HIST_VERBOSE);
    }
  }

  hist_save(hs, m->fnhis, HIST_ADDAHALF);
  fprintf(stderr, "simulation ended, doing WHAM\n");

  for ( iT = 0; iT < m->nT; iT++ ) {
    lj_close( lj[iT] );
  }

  free(mctot);
  free(mcacc);
  free(mcamp);
  return hs;
}



static hist_t *lj_domd(model_t *m, const double *beta)
{
  double emin, emax;
  int istep, iT;
  lj_t **lj, *ljtmp;
  hist_t *hs;

  /* round emin and emax to multiples of m->de */
  emin = (int) (m->nn * m->emin / m->de) * m->de;
  emax = (int) (m->nn * m->emax / m->de) * m->de;

  hs = hist_open(m->nT, emin, emax, m->de);

  /* randomize the initial state */
  mtscramble( time(NULL) );

  xnew(lj, m->nT);
  for ( iT = 0; iT < m->nT; iT++ ) {
    lj[iT] = lj_open(m->nn, m->rho, m->rcdef);
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
        for ( i = 0; i < m->nn; i++ ) {
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
  return hs;
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
    if ( m->simul == SIMUL_MD ) {
      hs = lj_domd(m, beta);
    } else {
      hs = lj_domc(m, beta);
    }
  }

  whamx(hs, beta, lnz,
      m->mdiis_nbases, m->mdiis_damp,
      m->mdiis_update_method, m->mdiis_threshold,
      m->itmax, m->tol, m->itmin, m->verbose,
      m->fnlndos, m->fneav, m->wham_method);

  if ( m->verbose ) {
    double lnzref0 = 0;

    for ( iT = 0; iT < hs->rows; iT++ ) {
      double tp = 1/beta[iT], tot, eav, var;
      double eref, pref, lnzref, muref;

      eav = hist_getave(hs, iT, &tot, &var);
      eref = ljeos3d_get(m->rho, tp, &pref, &lnzref, &muref);
      lnzref *= -m->nn / tp;
      if ( iT == 0 ) lnzref0 = lnzref;
      printf("%3d %10.7f %14.7f %8.0f %15.7f %14.7f | %10.3f %10.3f\n",
          iT, tp, lnz[iT], tot, eav, sqrt(var),
          lnzref - lnzref0, eref * m->nn);
    }
  }

  hist_close(hs);
  free(beta);
  free(lnz);
  return 0;
}

