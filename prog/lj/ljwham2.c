/* test program for the two-dimensional WHAM */
#define WHAM2_MDIIS
#include "../wham2.h"
#include "lj.h"


int n = 108;
double rho = 0.7;
double rcdef = 2.5;
double dt = 0.002; /* MD time step */
double thdt = 0.02; /* thermostat time step */
double pdt = 1e-5; /* barostat time step */
int nequil = 1000;
int nsteps = 10000;

int nt = 10;
int np = 10;
int ntp = 100;
double emin = -800, emax = 200, de = 1.0;
double vmin =  100, vmax = 400, dv = 1.0;
int itmax = 100000;
double tol = 1e-10;
int nbases = 5;
int verbose = 0;

const char *fnlndos = "lndos2.dat";
const char *fneav = "eav2.dat";
const char *fnhist = "hist2.dat";

enum { METHOD_DIRECT = 0, METHOD_MDIIS = 1 };
int method = METHOD_DIRECT;



int main(int argc, char **argv)
{
  hist2_t *hs;
  double *beta, *pres, *bp, *lnz, *epot, *vols, T, p;
  lj_t **lj;
  int it, ip, itp, istep;

  if ( argc > 1 ) { /* choose the method of WHAM */
    method = atoi(argv[1]);
  }

  ntp = nt * np;
  xnew(lj, ntp);
  xnew(beta, ntp);
  xnew(pres, ntp);
  xnew(bp, ntp);
  xnew(lnz, ntp);
  xnew(epot, ntp);
  xnew(vols, ntp);

  /* initialize the systems */
  itp = 0;
  for ( it = 0; it < nt; it++ ) {
    T = 1.0 + 0.1 * it;
    for ( ip = 0; ip < np; ip++ ) {
      p = 1.0 + 0.1 * ip;
      lj[itp] = lj_open(n, rho, rcdef);
      beta[itp] = 1./T;
      pres[itp] = p;
      bp[itp] = beta[itp] * pres[itp];
      lnz[itp] = epot[itp] = 0;
      itp++;
    }
  }
  hs = hist2_open(ntp, emin, emax, de, vmin, vmax, dv);

  /* try to load the histogram, if it fails, do simulations */
  if ( 0 != hist2_load(hs, fnhist, HIST2_VERBOSE) ) {
    /* do simulations */
    for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
      for ( itp = 0; itp < ntp; itp++ ) { /* loop over replicas */
        lj_vv(lj[itp], dt);
        lj[itp]->ekin = lj_vrescale(lj[itp], 1/beta[itp], thdt);
        if ( istep % 5 == 0 ) {
          lj_langp0(lj[itp], pdt, 1/beta[itp], pres[itp], 0);
        }
        epot[itp] = lj[itp]->epot;
        vols[itp] = lj[itp]->vol;
      }
      if ( istep % 1000 == 0 ) printf("t %d\n", istep);
      if ( istep <= nequil ) continue;
      hist2_add(hs, epot, vols, 1, 1.0, HIST2_VERBOSE);
    }

    hist2_save(hs, fnhist, HIST2_ADDAHALF);
    fprintf(stderr, "simulation ended, doing WHAM\n");
  }

  /* do WHAM */
  if ( method == METHOD_DIRECT ) {
    wham2(hs, beta, bp, lnz,
        itmax, tol, verbose, fnlndos, fneav);
  } else {
    wham2_mdiis(hs, beta, bp, lnz, nbases, 1.0,
        itmax, tol, verbose, fnlndos, fneav);
  }

  /* clean up */
  hist2_close(hs);
  for ( itp = 0; itp < ntp; itp++ )
    lj_close( lj[itp] );
  free(beta);
  free(pres);
  free(bp);
  free(lnz);
  free(epot);
  free(vols);
  return 0;
}

