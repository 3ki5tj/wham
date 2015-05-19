#ifndef MODEL_H__
#define MODEL_H__





#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>





enum { WHAM_DIRECT = 0, WHAM_MDIIS = 1, WHAM_NMETHODS };
const char *wham_methods[] = {"direct", "mdiis"};



typedef struct {
  char *prog;
  char *fnhis;
  double de;
  int wham_method;
  int itmax;
  double tol;
  int mdiis_nbases; /* number of bases in MDIIS */
  double mdiis_damp; /* mixing factor in MDIIS */
  char *fnlndos;
  char *fneav;
  double actmax; /* cutoff time of the autocorrelation function */
  double acmin; /* minimum of the autocorrelation function */
  char *fnac;
  char *fninp;
  int verbose;
} model_t;



/* set default values of the parameters */
__inline static void model_default(model_t *m)
{
  memset(m, 0, sizeof(*m));
  m->prog = NULL;
  m->fnhis = "hist.dat";
  m->de = 1.0;
  m->wham_method = WHAM_DIRECT;
  m->itmax = 100000;
  m->tol = 1e-10;
  m->mdiis_nbases = 5;
  m->mdiis_damp = 1.0;
  m->fnlndos = NULL;
  m->fneav = NULL;
  m->actmax = 0;
  m->acmin = 0;
  m->fninp = NULL;
  m->fnac = NULL;
  m->verbose = 0;
}



/* return the index of string from a predefined array */
__inline static int model_select(const char *s, int n, const char **arr)
{
  int i;

  for ( i = 0; i < n; i++ )
    if ( strcmp(arr[i], s) == 0 )
      return i;

  if ( isdigit(s[0]) ) {
    i = atoi(s);
    if ( i >= 0 && i < n ) {
      return i;
    }
  }

  fprintf(stderr, "Error: cannot select %s from the array of %d items:", s, n);
  for ( i = 0; i < n; i++ )
    fprintf(stderr, " %s", arr[i]);
  fprintf(stderr, "\n");

  exit(1);
  return 0;
}



/* print help message and die */
__inline static void model_help(const model_t *m)
{
  fprintf(stderr, "WHAM for xvg file\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [Options] input\n\n", m->prog);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -i:            set the input file, default: %s\n", m->fninp);
  fprintf(stderr, "  --fnhis=:      set the histogram file, default: %s\n", m->fnhis);
  fprintf(stderr, "  --de=:         set the energy bin size, default %g\n", m->de);
  fprintf(stderr, "  --wham=:       set the WHAM method, 'direct' or 'mdiis', default: %s\n", wham_methods[m->wham_method]);
  fprintf(stderr, "  --itmax=:      set the maximal number of iterations, default %d\n", m->itmax);
  fprintf(stderr, "  --tol=:        set the tolerance of error, default %g\n", m->tol);
  fprintf(stderr, "  --nbases=:     set the number of bases in the MDIIS method, default: %d\n", m->mdiis_nbases);
  fprintf(stderr, "  --mdamp=:      set the mixing factor in the MDIIS method, default: %g\n", m->mdiis_damp);
  fprintf(stderr, "  --fndos=:      set the file for the density of states, default %s\n", m->fnlndos);
  fprintf(stderr, "  --fneav=:      set the file for the average energy, default %s\n", m->fneav);
  fprintf(stderr, "  --actmax=:     set the time cutoff of the autocorrelation function, default: %g\n", m->actmax);
  fprintf(stderr, "  --acmin=:      set thet minimal value of the autocorrelation function, default: %g\n", m->acmin);
  fprintf(stderr, "  --fnac=:       set the autocorrelation functions, default: %s\n", m->fnac);
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc., default %d\n", m->verbose);
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
}



/* handle command line arguments */
__inline static void model_doargs(model_t *m, int argc, char **argv)
{
  int i, j, ch;
  char *p, *q;

  /* reset */
  model_default(m);

  /* set the program name */
  m->prog = argv[0];

  for ( i = 1; i < argc; i++ ) {
    /* it's an argument */
    if ( argv[i][0] != '-' ) {
      m->fninp = argv[i];
      continue;
    }

    /* long argument, like --help */
    if ( argv[i][1] == '-' ) {
      /* try to parse the argment
         e.g., `--prog=aaa' is parsed to `--prog' and `aaa' */
      p = argv[i] + 2;
      /* let q point to the argument of the option */
      if ( (q = strchr(p, '=')) != NULL ) {
        *q++ = '\0';
      } else {
        q = NULL;
      }

      if ( strcmp(p, "fnhis") == 0 ) {
        m->fnhis = q;
      } else if ( strcmp(p, "de") == 0 ) {
        m->de = atof(q);
      } else if ( strcmp(p, "wham") == 0 ) {
        m->wham_method = model_select(q, WHAM_NMETHODS, wham_methods);
      } else if ( strcmp(p, "itmax") == 0 ) {
        m->itmax = atoi(q);
      } else if ( strcmp(p, "tol") == 0 ) {
        m->tol = atof(q);
      } else if ( strcmp(p, "nbases") == 0 ) {
        m->mdiis_nbases = atoi(q);
      } else if ( strcmp(p, "mdamp") == 0 ) {
        m->mdiis_damp = atof(q);
      } else if ( strcmp(p, "fndos") == 0 ) {
        m->fnlndos = q;
      } else if ( strcmp(p, "fneav") == 0 ) {
        m->fneav = q;
      } else if ( strcmp(p, "actmax") == 0 ) {
        m->actmax = atof(q);
      } else if ( strcmp(p, "acmin") == 0 ) {
        m->acmin = atof(q);
      } else if ( strcmp(p, "fnac") == 0 ) {
        m->fnac = q;
      } else if ( strcmp(p, "help") == 0 ) {
        model_help(m);
      } else {
        fprintf(stderr, "unknown option %s\n", argv[i]);
        model_help(m);
      }

      continue;
    }

    /* it is an option
     * loop over characters in the options
     * in this way, `-vo' is understood as `-v -o' */
    for ( j = 1; (ch = argv[i][j]) != '\0'; j++ ) {
      if ( strchr("i", ch) != NULL ) {
        /* handle options that require an argument */
        q = p = argv[i] + j + 1;
        if ( *p != '\0' ) {
          /* the argument follows immediately after the option
           * e.g., -oa.dat */
          q = p;
        } else if ( ++i < argc ) {
          /* the option and argument are separated by a space
           * then the argument belongs to the next argv[] element,
           * hence ++i
           * e.g., -o a.dat */
          q = argv[i];
        } else {
          fprintf(stderr, "-%c requires an argument!\n", ch);
          model_help(m);
        }

        if ( ch == 'i' ) {
          m->fninp = q;
        }
        break; /* skip the rest of the characters in the option */
      } else if ( ch == 'v' ) {
        m->verbose++;
      } else if ( ch == 'h' ) {
        model_help(m);
      } else {
        fprintf(stderr, "unknown option %s, j %d, ch %c\n", argv[i], j, ch);
        model_help(m);
      }
    }
  }
}



#endif /* MODEL_H__ */
