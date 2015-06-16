#ifndef LSUTIL_H__
#define LSUTIL_H__





/* utilities to handle list files */





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



/* load the list of the names of energy files */
static char **getls(const char *fn,
    int *nbeta, double **beta)
{
  FILE *fp;
  int i;
  double tp, boltz;
  char buf[4096] = "", *p;
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
  fprintf(stderr, "%d temperatures\n", *nbeta);

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

    /* if the file extension is .xvg, then use
     * the GROMACS units, otherwise use the reduced units */
    p = strrchr(fns[i], '.');
    if ( p != NULL && strcmp(p, ".xvg") == 0 ) {
      boltz = BOLTZ;
    } else {
      boltz = 1;
    }

    (*beta)[i] = 1 / (boltz * tp);
  }

  fclose(fp);

  return fns;
}





#endif /* LSUTIL_H__ */
