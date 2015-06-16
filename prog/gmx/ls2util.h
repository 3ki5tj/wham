#ifndef LS2UTIL_H__
#define LS2UTIL_H__





/* utilities to handle two-dimensional list files */




/* GROMACS constants
 * according to src/gromacs/legacyheaders/physics.h */
#define BOLTZ   (1.380658e-23*6.0221367e23/1e3)
#define PRESFAC (16.6054) /* bar / pressure unity = 1e27/(6.0221367e23*100)
                             1 bar = 100 kPa = 100 kJ/m^3 = 10^2*6.0221367e23/1e27 kJ/mol/nm^3 */




/* determine the temperature, pressure and file name
 * from a line of the input list file */
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



/* load the list of the names of the energy-volume files
 * return the list of file names, and the inverse temperature
 * beta and beta*p values are saved in
 * (*beta)[] and (*bpres)[], respectively */
static char **getls(const char *fn,
    int *nbp, double **beta, double **bpres)
{
  FILE *fp;
  int i;
  double temp, pres, boltz, presfac;
  char buf[4096], *p;
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
  fprintf(stderr, "%d temperatures\n", *nbp);

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
  
    /* determine the unit system
     * if the ending is .xvg, we use the GROMACS unit system
     * otherwise, we use the reduced units */
    p = strrchr(fn, '.');
    if ( p != NULL && strcmp(p, ".xvg" ) == 0 ) {
      boltz = BOLTZ;
      presfac = PRESFAC;
    } else {
      boltz = 1;
      presfac = 1;
    }

    (*beta)[i] = 1 / (boltz * temp);
    (*bpres)[i] = (*beta)[i] * pres / presfac;
  }

  fclose(fp);

  return fns;
}





#endif /* LS2UTIL_H__ */

