/* test program for WHAM */
#define WHAM_MDIIS
#include "../wham.h"
#define IS2_LB 6
#include "is2.h"


int nequil = 100000;
int nsteps = 10000000;
int verbose = 0;

int ntp = 80;
double xmin = -2*IS2_N - 4, xmax = 4, dx = 4;
int itmax = 100000;
double tol = 1e-7;
int nbases = 5;

const char *fnlndos = "lndos.dat";
const char *fneav = "eav.dat";
const char *fnhist = "hist.dat";

enum { METHOD_DIRECT = 0, METHOD_MDIIS = 1 };
int method = METHOD_DIRECT;




int main(int argc, char **argv)
{
  hist_t *hs;
  double *beta, *lnz, *epot;
  is2_t **is;
  int itp, istep;

  if ( argc > 1 ) {
    method = atoi(argv[1]);
  }

  xnew(is, ntp);
  xnew(beta, ntp);
  xnew(lnz, ntp);
  xnew(epot, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    is[itp] = is2_open(IS2_L);
    beta[itp] = 1./(1.5 + 0.02 * itp);
    IS2_SETPROBA(is[itp], beta[itp]);
    lnz[itp] = epot[itp] = 0;
  }
  hs = hist_open(ntp, xmin, xmax, dx);

  /* try to load the histogram, if it fails, do simulations */
  if ( 0 != hist_load(hs, fnhist, HIST_VERBOSE) ) {
    int id, h;

    /* do the simulations */
    for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
      for ( itp = 0; itp < ntp; itp++ ) {
        IS2_PICK(is[itp], id, h);
        if ( h < 0 || mtrand() <= is[itp]->uproba[h] ) {
          IS2_FLIP(is[itp], id, h);
        }
        epot[itp] = is[itp]->E;
      }
      if ( istep <= nequil ) continue;
      //for ( itp = 0; itp < ntp; itp++ ) printf("itp %d, ep %d\n", itp, is[itp]->E);
      //printf("hs->xmin %g\n", hs->xmin); getchar();
      hist_add(hs, epot, 1, 0);
    }

    hist_save(hs, fnhist, 0);
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
    is2_close( is[itp] );
  free(beta);
  free(lnz);
  free(epot);
  return 0;
}

