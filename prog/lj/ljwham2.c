/* test program for the two-dimensional WHAM */
#define ENABLE_MDIIS
#include "../wham2.h"
#include "lj.h"
#include "ljeos.h"
#include <time.h>
#define WHAM
#define LJ_MODEL
#include "../whammodel.h"




static void model_default_lj2(model_t *m)
{
  model_default(m);
  m->nn = 256;
  m->rho = 0.3;
  m->rcdef = 1e+30;
  m->mcamp = 0.5;
  m->mddt = 0.002;
  m->thdt = 0.02;
  m->pdt = 1e-5;
  m->de = 1.0;
  m->emin = -10.0;
  m->emax = -0.0;
  m->dv = 2.0;
  m->vmin = 0.0;
  m->vmax = 25.0;
  m->nT = 6;
  m->Tmin = 1.2;
  m->Tdel = 0.1;
  m->nP = 3;
  m->Pmin = 0.1;
  m->Pdel = 0.05;
  m->nstadj = 20000;
  m->nequil = 100000;
  m->nsteps = 2000000;
  m->simul = SIMUL_MC;
  m->vamp = 0.1;
  m->nstvmov = 10;
  m->defsetup = 0;
}



static hist2_t *lj_domc(model_t *m, int ntp, double *beta, double *bp)
{
  double emin, emax, vmin, vmax;
  int istep, itp, jtp, acc;
  double *mctot, *mcacc, *mcamp;
  double *vtot, *vacc, *vamp;
  double retot = DBL_MIN, reacc = 0.0;
  lj_t **lj, *ljtmp;
  hist2_t *hs;

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
    double rho = bp[itp] / beta[itp];
    lj[itp] = lj_open(m->nn, rho, m->rcdef);
    lj[itp]->dof = m->nn * D;
    lj_energy( lj[itp] );
  }

  xnew(mctot, ntp);
  xnew(mcacc, ntp);
  xnew(mcamp, ntp);
  xnew(vtot, ntp);
  xnew(vacc, ntp);
  xnew(vamp, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    mctot[itp] = DBL_MIN;
    mcacc[itp] = 0;
    mcamp[itp] = m->mcamp;
    vtot[itp] = DBL_MIN;
    vacc[itp] = 0;
    vamp[itp] = m->vamp;
  }

  /* do simulations */
  for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
    for ( itp = 0; itp < ntp; itp++ ) { /* loop over replicas */
      acc = lj_metro(lj[itp], mcamp[itp], beta[itp]);
      mcacc[itp] += acc;
      mctot[itp] += 1;
      if ( istep % m->nstvmov == 0 ) {
        acc = lj_mcvmov(lj[itp], vamp[itp],
            1/beta[itp], bp[itp]/beta[itp], 0);
        //if ( acc ) {
        //  printf("%d, epot %g(%g,%g,%g), vol %g -> ", itp, lj[itp]->epot, lj[itp]->ep6, lj[itp]->ep12, lj[itp]->vir, lj[itp]->vol);
        //  lj_energy(lj[itp]);
        //  printf("epot %g(%g,%g,%g), vol %g\n", lj[itp]->epot, lj[itp]->ep6, lj[itp]->ep12, lj[itp]->vir, lj[itp]->vol); getchar();
        //}
        vacc[itp] += acc;
        vtot[itp] += 1;
      }
    }

    if ( m->re ) {
      /* replica exchange: randomly swap configurations */
      double dbdE, r;

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
        /* swap the models */
        ljtmp = lj[itp], lj[itp] = lj[jtp], lj[jtp] = ljtmp;
      }

      if ( istep > m->nequil ) {
        reacc += acc;
        retot += 1.0;
      }
    }

    if ( istep <= m->nequil ) {
      /* adjust the MC move amplitude */
      if ( istep % m->nstadj == 0 ) {
        for ( itp = 0; itp < ntp; itp++ ) {
          fprintf(stderr, "T %g, mcamp %g, mcacc %g%%, vamp %g, vacc %g%%\n",
              1/beta[itp], mcamp[itp], 100*mcacc[itp]/mctot[itp],
              vamp[itp], 100*vacc[itp]/vtot[itp]);
          mcamp[itp] *= sqrt( mcacc[itp] / mctot[itp] / 0.5 );
          vamp[itp] *= sqrt( vacc[itp] / vtot[itp] / 0.5 );
        }
      }
      continue;
    }

    for ( itp = 0; itp < ntp; itp++ ) {
      hist2_add1(hs, itp, lj[itp]->epot, lj[itp]->vol,
          1.0, HIST2_VERBOSE);
    }
  }

  hist2_save(hs, m->fnhis2, HIST2_ADDAHALF | HIST2_NOZEROES | HIST2_VERBOSE);
  fprintf(stderr, "simulation ended, doing WHAM, reacc = %g%%\n",
      100*reacc/retot);

  for ( itp = 0; itp < ntp; itp++ ) {
    lj_close( lj[itp] );
  }

  free(mctot);
  free(mcacc);
  free(mcamp);
  free(vtot);
  free(vacc);
  free(vamp);
  return hs;
}



static hist2_t *lj_domd(model_t *m, int ntp, double *beta, double *bp)
{
  double emin, emax, vmin, vmax;
  int istep, itp, jtp, acc;
  double retot = DBL_MIN, reacc = 0.0;
  lj_t **lj, *ljtmp;
  hist2_t *hs;

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
    double rho = bp[itp] / beta[itp];
    lj[itp] = lj_open(m->nn, rho, m->rcdef);
  }

  /* do simulations */
  for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
    for ( itp = 0; itp < ntp; itp++ ) { /* loop over replicas */
      lj_vv(lj[itp], m->mddt);
      lj[itp]->ekin = lj_vrescale(lj[itp], 1/beta[itp], m->thdt);
      if ( istep % m->nstvmov == 0 ) {
        lj_langp0(lj[itp], m->pdt, 1/beta[itp], bp[itp]/beta[itp], 0);
      }
    }

    if ( m->re ) {
      /* replica exchange: randomly swap configurations */
      double dbdE, r;

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

  hist2_save(hs, m->fnhis2, HIST2_ADDAHALF | HIST2_NOZEROES | HIST2_VERBOSE);
  fprintf(stderr, "simulation ended, doing WHAM, reacc = %g%%\n",
      100*reacc/retot);

  for ( itp = 0; itp < ntp; itp++ ) {
    lj_close( lj[itp] );
  }

  return hs;
}



/* return the density at which pressure is `pres` */
static double refsolve(double rho, double tp, double pres)
{
  double drho, pref, fref, muref;
  double rhol = 0, rhoh = DBL_MAX;

  while ( 1 ) {
    ljeos3d_get(rho, tp, &pref, &fref, &muref);
    drho = (pres - pref) / tp;
    if ( fabs(drho) < 1e-10 ) {
      break;
    }
    if ( pref > pres && rho < rhoh ) /* update rhoh */
      rhoh = rho;
    if ( pref < pres && rho > rhol ) /* update rhol */
      rhol = rho;
    rho += drho;
    if ( rho > rhoh || rho < rhol ) {
      rho = (rhoh + rhol) / 2;
    }
  }
  return rho;
}



int main(int argc, char **argv)
{
  model_t m[1];
  hist2_t *hs;
  double *beta, *pres, *bp, *lnz;
  int it, ip, itp, ntp;
#define NTPDEF 24
  double tpdef[NTPDEF]   = {1.30, 1.30, 1.30, 1.30, 1.30, 1.30, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60};
  double presdef[NTPDEF] = {0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.09, 0.10, 0.12, 0.15, 0.18, 0.21, 0.10, 0.12, 0.17, 0.22, 0.27, 0.32, 0.11, 0.12, 0.20, 0.28, 0.36, 0.44};

  model_default_lj2(m);
  model_doargs(m, argc, argv);

  if ( m->defsetup ) {
    ntp = NTPDEF;
  } else {
    ntp = m->nT * m->nP;
  }

  xnew(beta, ntp);
  xnew(pres, ntp);
  xnew(bp, ntp);
  xnew(lnz, ntp);

  /* set up the temperatures and pressures */
  if ( m->defsetup ) {
    for ( itp = 0; itp < ntp; itp++ ) {
      beta[itp] = 1/tpdef[itp];
      pres[itp] = presdef[itp];
    }
  } else {
    itp = 0;
    for ( it = 0; it < m->nT; it++ ) {
      double T = m->Tmin + m->Tdel * it;
      for ( ip = 0; ip < m->nP; ip++ ) {
        double P = m->Pmin + m->Pdel * ip;
        beta[itp] = 1./T;
        pres[itp] = P;
        itp++;
      }
    }
  }

  for ( itp = 0; itp < ntp; itp++ ) {
    bp[itp] = beta[itp] * pres[itp];
    lnz[itp] = 0;
  }

  if ( m->loadprev ) {
    /* load from existing histogram */
    if ( (hs = hist2_initf(m->fnhis2, HIST2_INT | HIST2_VERBOSE)) == NULL ) {
      return -1;
    }
    if ( hs->rows != ntp ) {
      fprintf(stderr, "number of histograms mismatch %d vs %d\n", hs->rows, ntp);
      return -1;
    }
  } else {
    if ( m->simul == SIMUL_MD ) {
      hs = lj_domd(m, ntp, beta, bp);
    } else {
      hs = lj_domc(m, ntp, beta, bp);
    }
  }

  /* do WHAM */
  wham2x(hs, beta, bp, lnz,
      m->damp, m->mdiis_nbases,
      m->mdiis_update_method, m->mdiis_threshold,
      m->itmin, m->itmax, m->tol, m->verbose,
      m->fnlndos2, m->fneav2, m->wham_method);

  if ( m->verbose ) {
    double lnzref0 = 0;

    for ( itp = 0; itp < ntp; itp++ ) {
      double tp = 1/beta[itp], rho;
      double tot, eav, vav, see, sev, svv, sige, sigv;
      double eref, pref, fref, muref, lnzref;

      eav = hist2_getave(hs, itp, &tot, &vav, &see, &sev, &svv);
      rho = m->nn / vav;
      /* find the density at which P = pres[itp] */
      rho = refsolve(rho, tp, pres[itp]);
      eref = ljeos3d_get(rho, tp, &pref, &fref, &muref);
      //printf("tp %g, beta %g, bp %g, %g, -bfref %g\n", tp, beta[itp], beta[itp]*pres[itp], beta[itp]*pref, -beta[itp]*fref*m->nn);
      lnzref = -m->nn * fref / tp;
      if ( itp == 0 ) lnzref0 = lnzref;
      sige = sqrt(see);
      sigv = sqrt(svv);
      printf("%3d %10.7f %10.7f %14.7f %8.0f "
          "%15.7f %14.7f %15.7f %14.7f %15.7f | "
          "%10.3f %10.3f %10.3f\n",
          itp, tp, pres[itp], lnz[itp], tot,
          eav, sige, vav, sigv, sev,
          lnzref - lnzref0, eref * m->nn, m->nn / rho);
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

