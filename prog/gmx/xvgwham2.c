/* WHAM for a set of xvg files */
#include "xvgmodel.h"
#include "xvg.h"
#define WHAM2_MDIIS
#include "../wham2.h"





/* according to src/gromacs/legacyheaders/physics.h */
#define BOLTZ   (1.380658e-23*6.0221367e23/1e3)
#define PRESFAC (16.6054) /* bar / pressure unity = 1e27/(6.0221367e23*100) */




/* determine the temperature, pressure and file name
 * from the input line */
static double parsefn(char *buf, double *pres, char *fn)
{
  char *p, *q;
  double temp = 300;

  *pres = 1;

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
    p++;
    for ( q = p; *q != '\0' && *q != '\n'; q++ ) {
      if ( isspace(*q) ) {
        break;
      }
    }

    if ( *q != '\0' && *q != '\n' ) {
      sscanf(buf, "%lf %lf %s", &temp, pres, fn);
      return temp;
    }
  } else {
    p = buf;
  }

  /* copy the file name */
  strcpy(fn, p);

  /* try to determine the temperature and pressure from the file name */
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

  sscanf(p, "T%lfP%lf", &temp, pres);

  return temp;
}



/* load the list of energy files */
static char **getls(const char *fn,
    int *nbp, double **beta, double **bpres)
{
  FILE *fp;
  int i;
  double temp, pres;
  char buf[4096];
  char **fns;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return NULL;
  }

  /* determine the number of lines */
  *nbp = 0;
  while ( fgets(buf, sizeof buf, fp) ) {
    if ( buf[0] == '#' || isspace(buf[0]) ) { /* skip a comment */
      continue;
    }
    *nbp += 1;
  }
  if ( *nbp <= 1 ) {
    fprintf(stderr, "insufficient number of lines, %d\n", *nbp);
    return NULL;
  }
  printf("%d temperatures\n", *nbp);

  /* go back to the beginning of the file */
  rewind(fp);

  xnew(fns, *nbp);
  xnew(*beta, *nbp);
  xnew(*bpres, *nbp);
  for ( i = 0; i < *nbp; i++ ) {
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
    temp = parsefn(buf, &pres, fns[i]);
    (*beta)[i] = 1 / (BOLTZ * temp);
    (*bpres)[i] = (*beta)[i] * pres / PRESFAC;
  }

  fclose(fp);

  return fns;
}



/* construct the histogram */
static hist2_t *mkhist2(const char *fnls,
    double **beta, double **bpres,
    double de, double dv, const char *fnhis)
{
  hist2_t *hs;
  int i, j, nbp, k;
  char **fns;
  xvg_t **xvg = NULL;
  double evmin[3] = {1e30, 1e30, 1e30}, evmax[3] = {-1e30, -1e30, -1e30};
  double evmin1[3], evmax1[3];

  if ( (fns = getls(fnls, &nbp, beta, bpres)) == NULL ) {
    return NULL;
  }

  xnew(xvg, nbp);
  for ( i = 0; i < nbp; i++ ) {
    if ( (xvg[i] = xvg_load(fns[i])) == NULL ) {
      return NULL;
    }
    xvg_minmax(xvg[i], evmin1, evmax1);
    for ( k = 0; k < 3; k++ ) {
      if ( evmin1[k] < evmin[k] ) {
        evmin[k] = evmin1[k];
      }
      if ( evmax1[k] > evmax[k] ) {
        evmax[k] = evmax1[k];
      }
    }
  }

  evmin[0] = ((int) (evmin[0] / de) - 1) * de;
  evmax[0] = ((int) (evmax[0] / de) + 1) * de;
  evmin[2] = ((int) (evmin[2] / dv) - 1) * dv;
  evmax[2] = ((int) (evmax[2] / dv) + 1) * dv;

  hs = hist2_open(nbp, evmin[0], evmax[0], de,
      evmin[2], evmax[2], dv);

  for ( i = 0; i < nbp; i++ ) {
    for ( j = 0; j < xvg[i]->n; j++ ) {
      hist2_add1(hs, i,
          xvg[i]->y[0][j], xvg[i]->y[2][j], 1.0, HIST2_VERBOSE);
    }
  }

  hist2_save(hs, fnhis, HIST2_ADDAHALF|HIST2_VERBOSE);

  for ( i = 0; i < nbp; i++ ) {
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
  hist2_t *hs = NULL;
  int i, nbp;
  double *beta, *bpres, *lnz;

  model_default(m);
  m->de = 10.0;
  model_doargs(m, argc, argv);

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  /* try to load from the existing histogram */
  if ( m->loadprev ) {
    hs = hist2_initf(m->fnhis2);
    if ( hs == NULL ) {
      return -1;
    }

    getls(m->fninp, &nbp, &beta, &bpres);
    if ( nbp != hs->rows ) {
      fprintf(stderr, "%s: different rows %d vs %d\n",
          m->fnhis2, nbp, hs->rows);
      hist2_close(hs);
      return -1;
    }
  } else {
    hs = mkhist2(m->fninp, &beta, &bpres, m->de, m->dv, m->fnhis2);
    if ( hs == NULL ) {
      return -1;
    }
  }

  xnew(lnz, hs->rows);
  for ( i = 0; i < hs->rows; i++ ) {
    lnz[i] = 0;
  }

  if ( m->wham_method == WHAM_DIRECT ) {
    wham2(hs, beta, bpres, lnz,
        m->itmax, m->tol, m->verbose, m->fnlndos, m->fneav);
  } else {
    wham2_mdiis(hs, beta, bpres, lnz,
        m->mdiis_nbases, m->mdiis_damp,
        m->itmax, m->tol, m->verbose, m->fnlndos, m->fneav);
  }

  hist2_close(hs);
  free(beta);
  free(lnz);
  return 0;
}

