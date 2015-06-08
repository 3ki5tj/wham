/* MBAR for a set of xvg files */
#include <time.h>
#include "../mtrand.h"
#define MTRAND
#include "xvg.h"
#define MBAR_ENABLE_MDIIS
#include "mbar.h"
#include "../whammodel.h"





/* according to src/gromacs/legacyheaders/physics.h */
#define BOLTZ (1.380658e-23*6.0221367e23/1e3)





/* determine the temperature and file name
 * from a line of the input list file */
static double parsefn(char *buf, char *fn)
{
  char *p, *q;
  double temp = 300;

  /* remove trailing spaces */
  for ( p = buf + strlen(buf) - 1; isspace(*p); p-- ) {
    *p = '\0';
  }

  /* check if there is a space */
  for ( p = buf; *p != '\0' && *p != '\n'; p++ ) {
    if ( isspace(*p) ) {
      break;
    }
  }

  /* there is a space */
  if ( *p != '\0' && *p != '\n' ) {
    sscanf(buf, "%lf %s", &temp, fn);
    return temp;
  }

  strcpy(fn, buf);

  /* try to determine the temperature from the file name */
  q = strrchr(buf, '/');
  if ( q == NULL ) return temp;
  *q = '\0';

  /* let p point to the directory containing the file */
  p = strrchr(buf, '/');
  if ( p == NULL ) {
    p = buf;
  } else {
    p++;
  }

  /* skip the character 'T' */
  if ( *p == 'T' ) {
    p++;
  }
  sscanf(p, "%lf", &temp);

  return temp;
}



/* load the list of energy files */
static char **getls(const char *fn,
    int *nbeta, double **beta)
{
  FILE *fp;
  int i;
  double tp;
  char buf[4096] = "";
  char **fns;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return NULL;
  }

  /* determine the number of lines */
  *nbeta = 0;
  while ( fgets(buf, sizeof buf, fp) ) {
    if ( buf[0] == '#' || isspace(buf[0]) ) { /* skip a comment */
      continue;
    }
    *nbeta += 1;
  }
  if ( *nbeta <= 1 ) {
    fprintf(stderr, "insufficient number of lines, %d\n", *nbeta);
    return NULL;
  }
  printf("%d temperatures\n", *nbeta);

  /* go back to the beginning of the file */
  rewind(fp);

  xnew(fns,   *nbeta);
  xnew(*beta, *nbeta);
  for ( i = 0; i < *nbeta; i++ ) {
    while ( 1 ) { /* loop till we get a non-comment line */
      if ( fgets(buf, sizeof buf, fp) == NULL ) {
        fprintf(stderr, "cannot read line %d from %s\n", i, fn);
        fclose(fp);
        return NULL;
      }
      if ( buf[0] != '#' ) {
        break;
      }
    }

    /* copy the line */
    xnew(fns[i], strlen(buf) + 1);
    tp = parsefn(buf, fns[i]);
    (*beta)[i] = 1 / (BOLTZ * tp);
  }

  fclose(fp);

  return fns;
}



/* construct the histogram */
static xvg_t **mkxvg(const char *fnls,
    int *nbeta, double **beta, double de,
    double radd)
{
  int i;
  char **fns;
  xvg_t **xvg = NULL;
  double xmin[2], xmax[2], emin = 1e30, emax = -1e30;

  /* scramble the random number seed */
  mtscramble( time(NULL) );

  if ( (fns = getls(fnls, nbeta, beta)) == NULL ) {
    return NULL;
  }

  xnew(xvg, *nbeta);
  for ( i = 0; i < *nbeta; i++ ) {
    xvg[i] = xvg_load(fns[i], radd);
    xvg_minmax(xvg[i], xmin, xmax);
    if ( xmin[0] < emin ) {
      emin = xmin[0];
    }
    if ( xmax[0] > emax ) {
      emax = xmax[0];
    }
  }

  emin = ((int) (emin / de) - 1) * de;
  emax = ((int) (emax / de) + 1) * de;

  for ( i = 0; i < *nbeta; i++ ) {
    free(fns[i]);
  }
  free(fns);

  return xvg;
}



int main(int argc, char **argv)
{
  model_t m[1];
  xvg_t **xvg = NULL;
  int i, nbeta;
  double *beta, *lnz;

  model_default(m);
  model_doargs(m, argc, argv);

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  xvg = mkxvg(m->fninp, &nbeta, &beta, m->de, m->radd);
  if ( xvg == NULL ) {
    return -1;
  }

  xnew(lnz, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    lnz[i] = 0;
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    mbar(nbeta, xvg, beta, lnz,
        m->itmax, m->tol, m->verbose);
  } else {
    mbar_mdiis(nbeta, xvg, beta, lnz,
        m->mdiis_nbases, m->mdiis_damp,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmax, m->tol, m->verbose);
  }

  for ( i = 0; i < nbeta; i++ ) {
    fprintf(stderr, "%3d %8.5f %10.3f %8d\n",
        i, beta[i], lnz[i], xvg[i]->n);
  }

  for ( i = 0; i < nbeta; i++ ) {
    free(xvg[i]);
  }
  free(xvg);
  free(beta);
  free(lnz);
  return 0;
}

