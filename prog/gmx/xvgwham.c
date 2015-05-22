/* WHAM for a set of xvg files */
#include "../whammodel.h"
#include "xvg.h"
#define WHAM_MDIIS
#include "../wham.h"
#include <time.h>
#include "../mtrand.h"





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
  char buf[4096];
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
static hist_t *mkhist(const char *fnls,
    double **beta, double de,
    const char *fnhis, double radd)
{
  hist_t *hs;
  int i, j, nbeta;
  char **fns;
  xvg_t **xvg = NULL;
  double emin = 1e30, emax = -1e30, emin1 = 1e30, emax1 = -1e30;

  /* scramble the random number seed */
  mtscramble( time(NULL) );

  if ( (fns = getls(fnls, &nbeta, beta)) == NULL ) {
    return NULL;
  }

  xnew(xvg, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    xvg[i] = xvg_load(fns[i]);
    xvg_minmax(xvg[i], &emin1, &emax1);
    if ( emin1 < emin ) {
      emin = emin1;
    }
    if ( emax1 > emax ) {
      emax = emax1;
    }
  }

  emin = ((int) (emin / de) - 1) * de;
  emax = ((int) (emax / de) + 1) * de;

  hs = hist_open(nbeta, emin, emax, de);

  for ( i = 0; i < nbeta; i++ ) {
    for ( j = 0; j < xvg[i]->n; j++ ) {
      if ( radd >= 1.0 || rand01() < radd ) {
        hist_add1(hs, i, xvg[i]->y[0][j], 1.0, HIST_VERBOSE);
      }
    }
  }

  hist_save(hs, fnhis, HIST_ADDAHALF|HIST_VERBOSE);

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
  hist_t *hs = NULL;
  int i, nbeta;
  double *beta, *lnz;

  model_default(m);
  model_doargs(m, argc, argv);

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  if ( m->loadprev ) {
    /* load from the existing histogram */
    if ( (hs = hist_initf(m->fnhis)) == NULL ) {
      return -1;
    }

    getls(m->fninp, &nbeta, &beta);
    if ( nbeta != hs->rows ) {
      fprintf(stderr, "%s: different rows %d vs %d\n",
          m->fnhis, nbeta, hs->rows);
      hist_close(hs);
      return -1;
    }
  } else {
    hs = mkhist(m->fninp, &beta, m->de, m->fnhis, m->radd);
    if ( hs == NULL ) {
      return -1;
    }
  }

  xnew(lnz, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) {
    lnz[i] = 0;
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

