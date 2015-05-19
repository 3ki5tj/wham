/* WHAM for a set of xvg files */
#include "model.h"
#include "xvg.h"
#define WHAM_MDIIS
#include "../wham.h"





/* according src/gromacs/legacyheaders/physics.h */
#define BOLTZ (1.380658e-23*6.0221367e23/1e3)





/* load the file list */
static char **getls(const char *fn,
    int *nbeta, double **beta)
{
  FILE *fp;
  int i;
  double tp;
  char buf[4096];
  char **fns;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return NULL;
  }

  if ( fgets(buf, sizeof buf, fp) == NULL
    || buf[0] != '#' ) {
    fprintf(stderr, "no leading line in %s\n", fn);
    fclose(fp);
    return NULL;
  }

  sscanf(buf+1, "%d", nbeta);

  xnew(fns,   *nbeta);
  xnew(*beta, *nbeta);

  for ( i = 0; i < *nbeta; i++ ) {
    if ( fgets(buf, sizeof buf, fp) == NULL ) {
      fprintf(stderr, "cannot read line %d from %s\n", i, fn);
      fclose(fp);
      return NULL;
    }

    /* copy the line */
    xnew(fns[i], strlen(buf) + 1);
    sscanf(buf, "%lf %s", &tp, fns[i]);
    (*beta)[i] = 1 / (BOLTZ * tp);
  }

  fclose(fp);

  return fns;
}



/* construct the histogram */
static hist_t *mkhist(model_t *m, double **beta)
{
  hist_t *hs;
  int i, j, nbeta;
  char **fns;
  xvg_t **xvg;
  double emin = 1e30, emax = -1e30, emin1 = 1e30, emax1 = -1e30;

  fns = getls(m->fninp, &nbeta, beta);
  xnew(xvg, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    xvg[i] = xvg_load(fns[i]);
    emin1 = xvg_minmax(xvg[i], &emax1);
    if ( emin1 < emin ) {
      emin = emin1;
    }
    if ( emax1 > emax ) {
      emax = emax1;
    }
  }

  emin = ((int) (emin / m->de) - 1) * m->de;
  emax = ((int) (emax / m->de) + 1) * m->de;

  hs = hist_open(nbeta, emin, emax, m->de);

  for ( i = 0; i < nbeta; i++ ) {
    for ( j = 0; j < xvg[i]->n; j++ ) {
      hist_add1(hs, i, xvg[i]->y[j], 1.0, HIST_VERBOSE);
    }
  }

  hist_save(hs, m->fnhis, HIST_ADDAHALF|HIST_VERBOSE);

  for ( i = 0; i < nbeta; i++ ) {
    xvg_close(xvg[i]);
    free(fns[i]);
  }
  free(xvg);
  free(fns);

  return hs;
}



int main(int argc, char **argv)
{
  model_t m[1];
  hist_t *hs;
  int i;
  double *beta, *lnz;

  model_default(m);
  model_doargs(m, argc, argv);

  hs = mkhist(m, &beta);
  xnew(lnz, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) {
    lnz[i] = 0;
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    wham(hs, beta, lnz,
        m->itmax, m->tol, m->fnlndos, m->fneav);
  } else {
    wham_mdiis(hs, beta, lnz, m->mdiis_nbases, 1.0,
        m->itmax, m->tol, 0, m->fnlndos, m->fneav);
  }

  hist_close(hs);
  free(beta);
  free(lnz);
  return 0;
}

