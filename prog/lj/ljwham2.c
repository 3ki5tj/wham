/* test program for the two-dimensional WHAM */
#define ENABLE_MDIIS
#include "../wham2.h"
#include "lj.h"
#include <time.h>
#define WHAM
#define LJ_MODEL
#include "../whammodel.h"




static void model_default_lj2(model_t *m)
{
  model_default(m);
  m->nn = 108;
  m->rho = 0.3;
  m->rcdef = 1e+30;
  m->mddt = 0.002;
  m->thdt = 0.02;
  m->pdt = 1e-5;
  m->de = 1.0;
  m->emin = -10.0;
  m->emax = -0.0;
  m->dv = 0.5;
  m->vmin = 0.0;
  m->vmax = 25.0;
  m->nT = 5;
  m->Tmin = 1.0;
  m->Tdel = 0.1;
  m->nP = 10;
  m->Pmin = 0.2;
  m->Pdel = 0.2;
  m->nequil = 10000;
  m->nsteps = 100000;
}



hist2_t *lj_simul(model_t *m, int ntp, double *beta, double *bp)
{
  hist2_t *hs;
  double emin, emax, vmin, vmax;
  int istep, itp;
  double retot = DBL_MIN, reacc = 0.0;
  lj_t **lj;

  /* round emin and emax to multiples of m->de */
  emin = (int) (m->nn * m->emin / m->de) * m->de;
  emax = (int) (m->nn * m->emax / m->de) * m->de;

  /* round vmin and vmax to multiples of m->dv */
  vmin = (int) (m->nn * m->vmin / m->dv) * m->dv;
  vmax = (int) (m->nn * m->vmax / m->dv) * m->dv;

  hs = hist2_open(ntp, emin, emax, m->de, vmin, vmax, m->dv);

  /* randomize the initial state */
  mtscramble( time(NULL) );

  xnew(lj, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    lj[itp] = lj_open(m->nn, m->rho, m->rcdef);
  }

  /* do simulations */
  for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
    for ( itp = 0; itp < ntp; itp++ ) { /* loop over replicas */
      lj_vv(lj[itp], m->mddt);
      lj[itp]->ekin = lj_vrescale(lj[itp], 1/beta[itp], m->thdt);
      if ( istep % 5 == 0 ) {
        lj_langp0(lj[itp], m->pdt, 1/beta[itp], bp[itp]/beta[itp], 0);
      }
    }

    if ( m->re ) {
      /* replica exchange: randomly swap configurations */
      int jtp, acc;
      double dbdE, r;
      lj_t *ljtmp;

      itp = (int) (rand01() * ntp);
      jtp = (itp + 1 + (int) (rand01() * (ntp - 1))) % ntp;
      dbdE = (beta[itp] - beta[jtp]) * (lj[itp]->epot - lj[jtp]->epot)
           + (bp[itp] - bp[jtp])* (lj[itp]->vol - lj[jtp]->vol);

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
        double scl = sqrt( beta[itp]/beta[jtp] );
        int i;

        /* scale the velocities */
        for ( i = 0; i < m->nn; i++ ) {
          vsmul(lj[itp]->v[i], scl);
          vsmul(lj[jtp]->v[i], 1/scl);
        }
        /* swap the models */
        ljtmp = lj[itp], lj[itp] = lj[jtp], lj[jtp] = ljtmp;
      }

      if ( istep > m->nequil ) {
        reacc += acc;
        retot += 1.0;
      }
    }

    if ( istep <= m->nequil ) continue;

    for ( itp = 0; itp < ntp; itp++ ) {
      hist2_add1(hs, itp, lj[itp]->epot, lj[itp]->vol,
          1.0, HIST2_VERBOSE);
    }
  }

  hist2_save(hs, m->fnhis2, HIST2_ADDAHALF);
  fprintf(stderr, "simulation ended, doing WHAM, reacc = %g%%\n",
      100*reacc/retot);

  for ( itp = 0; itp < ntp; itp++ ) {
    lj_close( lj[itp] );
  }

  return hs;
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
    if ( (hs = hist2_initf(m->fnhis2, HIST2_INT | HIST2_VERBOSE)) == NULL ) {
      return -1;
    }
  } else {
    hs = lj_simul(m, ntp, beta, bp);
  }

  /* do WHAM */
  wham2x(hs, beta, bp, lnz,
      m->mdiis_nbases, m->mdiis_damp,
      m->mdiis_update_method, m->mdiis_threshold,
      m->itmax, m->tol, m->itmin, m->verbose,
      m->fnlndos2, m->fneav2, m->wham_method);

  if ( m->verbose ) {
    for ( itp = 0; itp < hs->rows; itp++ ) {
      double tot, eav, vav, see, sev, svv;
      eav = hist2_getave(hs, itp, &tot, &vav, &see, &sev, &svv);
      printf("%3d %10.7f %10.7f %14.7f %8.0f "
          "%15.7f %14.7f %15.7f %14.7f %15.7f\n",
          itp, beta[itp], pres[itp], lnz[itp], tot,
          eav, sqrt(see), vav, sqrt(svv), sev);
    }
  }

  /* clean up */
  hist2_close(hs);
  free(beta);
  free(pres);
  free(bp);
  free(lnz);
  return 0;
}

