/* test program for WHAM */
#define WHAM_MDIIS
#include "../wham.h"
#include "lj.h"


int n = 108;
double rho = 0.3;
double rcdef = 2.5;
double dt = 0.002; /* MD time step */
double thdt = 0.02; /* thermostat time step */
int nequil = 4000;
int nsteps = 40000;

int ntp = 5;
double xmin = -600, xmax = 100, dx = 0.1;
int itmax = 100000;
double tol = 1e-10;
int nbases = 5;
int verbose = 0;

const char *fnlndos = "lndos.dat";
const char *fneav = "eav.dat";
const char *fnhist = "hist.dat";

enum { METHOD_DIRECT = 0, METHOD_MDIIS = 1 };
int method = METHOD_DIRECT;



int main(int argc, char **argv)
{
  hist_t *hs;
  double *beta, *lnz, *epot;
  lj_t **lj;
  int itp, istep;

  if ( argc > 1 ) {
    method = atoi(argv[1]);
  }

  xnew(lj, ntp);
  xnew(beta, ntp);
  xnew(lnz, ntp);
  xnew(epot, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    lj[itp] = lj_open(n, rho, rcdef);
    beta[itp] = 1./(0.8 + .3 * itp);
    lnz[itp] = epot[itp] = 0;
  }
  hs = hist_open(ntp, xmin, xmax, dx);

  /* try to load the histogram, if it fails, do simulations */
  if ( 0 != hist_load(hs, fnhist, HIST_VERBOSE) ) {
    /* do the simulations */
    for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
      for ( itp = 0; itp < ntp; itp++ ) {
        lj_vv(lj[itp], dt);
        lj[itp]->ekin = lj_vrescale(lj[itp], 1/beta[itp], thdt);
        epot[itp] = lj[itp]->epot;
      }
      if ( istep <= nequil ) continue;
      hist_add(hs, epot, 1, 0);
    }

    hist_save(hs, fnhist, HIST_ADDAHALF);
    fprintf(stderr, "simulation ended, doing WHAM\n");
  }

  if ( method == METHOD_DIRECT ) {
    wham(hs, beta, lnz,
        itmax, tol, verbose, fnlndos, fneav);
  } else {
    wham_mdiis(hs, beta, lnz, nbases, 1.0,
        itmax, tol, verbose, fnlndos, fneav);
  }
  hist_close(hs);
  for ( itp = 0; itp < ntp; itp++ )
    lj_close( lj[itp] );
  free(beta);
  free(lnz);
  free(epot);
  return 0;
}

