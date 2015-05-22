/* WHAM for two-dimensional Ising model */
#define IS2_MODEL
#include "../whammodel.h"
#define WHAM_MDIIS
#include "../wham.h"
#define IS2_LB 6
#include "is2.h"
#include <time.h>



double xmin = -2*IS2_N - 2, xmax = 2;



int main(int argc, char **argv)
{
  model_t m[1];
  hist_t *hs;
  double *beta, *lnz;
  int itp;

  model_default(m);
  m->de = 4;
  model_doargs(m, argc, argv);

  xnew(beta, m->ntp);
  xnew(lnz, m->ntp);
  for ( itp = 0; itp < m->ntp; itp++ ) {
    beta[itp] = 1./(m->tpmin + m->tpdel * itp);
    lnz[itp] = 0;
  }

  if ( m->loadprev ) {
    /* load from existing histogram */
    if ( (hs = hist_initf(m->fnhis)) == NULL ) {
      return -1;
    }
  } else {
    int id, h, istep;
    double *epot;
    is2_t **is;

    hs = hist_open(m->ntp, xmin, xmax, m->de);

    xnew(is, m->ntp);
    xnew(epot, m->ntp);
    for ( itp = 0; itp < m->ntp; itp++ ) {
      is[itp] = is2_open(IS2_L);
      IS2_SETPROBA(is[itp], beta[itp]);
      epot[itp] = 0;
    }

    /* randomize the initial state */
    mtscramble( time(NULL) );

    /* do the simulations */
    for ( istep = 1; istep <= m->nequil + m->nsteps; istep++ ) {
      /* MC for each temperature */
      for ( itp = 0; itp < m->ntp; itp++ ) {
        IS2_PICK(is[itp], id, h);
        if ( h < 0 || mtrand() <= is[itp]->uproba[h] ) {
          IS2_FLIP(is[itp], id, h);
        }
        epot[itp] = is[itp]->E;
      }

      if ( m->re ) {
        /* replica exchange: randomly swap configurations of
         * two neighboring temperatures */
        int jtp, acc;
        double dbdE, r;
        is2_t *istmp;
        unsigned utmp;

        itp = (int) (rand01() * (m->ntp - 1));
        jtp = itp + 1;
        dbdE = (beta[itp] - beta[jtp]) * (epot[itp] - epot[jtp]);
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
          istmp = is[itp], is[itp] = is[jtp], is[jtp] = istmp;
          /* swap the transition probabilities */
          utmp = is[itp]->uproba[2], is[itp]->uproba[2] = is[jtp]->uproba[2], is[jtp]->uproba[2] = utmp;
          utmp = is[itp]->uproba[4], is[itp]->uproba[4] = is[jtp]->uproba[4], is[jtp]->uproba[4] = utmp;
          /* swap the potential energies */
          r = epot[itp], epot[itp] = epot[jtp], epot[jtp] = r;
        }
      }

      if ( istep <= m->nequil ) continue;
      //for ( itp = 0; itp < m->ntp; itp++ ) printf("itp %d, ep %d\n", itp, is[itp]->E);
      //printf("hs->xmin %g\n", hs->xmin); getchar();
      hist_add(hs, epot, 1, 0);
    }

    hist_save(hs, m->fnhis, HIST_ADDAHALF);
    fprintf(stderr, "simulation ended %d steps, doing WHAM\n", m->nsteps);

    for ( itp = 0; itp < m->ntp; itp++ ) {
      is2_close( is[itp] );
    }
    free(epot);
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

