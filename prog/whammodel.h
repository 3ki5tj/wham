#ifndef WHAMMODEL_H__
#define WHAMMODEL_H__





#include "util.h"





enum { WHAM_DIRECT = 0, WHAM_MDIIS = 1, WHAM_NMETHODS };
const char *wham_methods[] = {"Direct", "MDIIS"};



typedef struct {
  char *prog;
  char *fnhis;
  char *fnhis2;
  int loadprev; /* load previous histogram */
  double radd; /* rate of adding trajectory frames into the histogram */
  double de;
  double dv;
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
  int re;
#ifdef IS2_MODEL
  int ntp;
  double tpmin;
  double tpdel;
  int nequil;
  int nsteps;
#endif
} model_t;



/* set default values of the parameters */
__inline static void model_default(model_t *m)
{
  memset(m, 0, sizeof(*m));
  m->prog = NULL;
  m->fnhis = "hist.dat";
  m->fnhis2 = "hist2.dat";
  m->loadprev = 0;
  m->radd = 1.0;
  m->de = 1.0;
  m->dv = 0.02;
  m->wham_method = WHAM_DIRECT;
  m->itmax = 100000;
  m->tol = 1e-10;
  m->mdiis_nbases = 10;
  m->mdiis_damp = 1.0;
  m->fnlndos = NULL;
  m->fneav = NULL;
  m->actmax = 0;
  m->acmin = 0;
  m->fninp = NULL;
  m->fnac = NULL;
  m->verbose = 0;
  m->re = 0;
#ifdef IS2_MODEL
  m->ntp = 80;
  m->tpmin = 1.5;
  m->tpdel = 0.02;
  m->nequil = 100000;
  m->nsteps = 10000000;
#endif
}



/* return the index of string from a predefined array */
__inline static int model_select(const char *s, int n, const char **arr)
{
  int i;

  for ( i = 0; i < n; i++ ) {
    if ( strcmpfuzzy(arr[i], s) == 0 ) {
      return i;
    }
  }

  if ( isdigit(s[0]) ) {
    i = atoi(s);
    if ( i >= 0 && i < n ) {
      return i;
    }
  }

  fprintf(stderr, "Error: cannot select %s from the array of %d items:", s, n);
  for ( i = 0; i < n; i++ ) {
    fprintf(stderr, " %s", arr[i]);
  }
  fprintf(stderr, "\n");

  exit(1);
  return 0;
}



/* print help message and die */
__inline static void model_help(const model_t *m)
{
  fprintf(stderr, "Weighted histogram analysis method (WHAM)\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [Options] input\n\n", m->prog);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -i:            set the input file, default: %s\n", m->fninp);
  fprintf(stderr, "  --fnhis=:      set the histogram file, default: %s\n", m->fnhis);
  fprintf(stderr, "  --fnhis2=:     set the 2D histogram file, default: %s\n", m->fnhis2);
  fprintf(stderr, "  -H:            load the previous histogram\n");
  fprintf(stderr, "  -r:            set the rate of adding trajectory frames into the histogram\n");
  fprintf(stderr, "  --de=:         set the energy bin size, default %g\n", m->de);
  fprintf(stderr, "  --dv=:         set the volume bin size, default %g\n", m->dv);
  fprintf(stderr, "  --wham=:       set the WHAM method, 'Direct' or 'MDIIS', default: %s\n", wham_methods[m->wham_method]);
  fprintf(stderr, "  --itmax=:      set the maximal number of iterations, default %d\n", m->itmax);
  fprintf(stderr, "  --tol=:        set the tolerance of error, default %g\n", m->tol);
  fprintf(stderr, "  --nbases=:     set the number of bases in the MDIIS method, default: %d\n", m->mdiis_nbases);
  fprintf(stderr, "  --mdamp=:      set the mixing factor in the MDIIS method, default: %g\n", m->mdiis_damp);
  fprintf(stderr, "  --fndos=:      set the file for the density of states, default %s\n", m->fnlndos);
  fprintf(stderr, "  --fneav=:      set the file for the average energy, default %s\n", m->fneav);
  fprintf(stderr, "  --actmax=:     set the time cutoff of the autocorrelation function, default: %g\n", m->actmax);
  fprintf(stderr, "  --acmin=:      set thet minimal value of the autocorrelation function, default: %g\n", m->acmin);
  fprintf(stderr, "  --fnac=:       set the autocorrelation functions, default: %s\n", m->fnac);
  fprintf(stderr, "  --re:          do replica exchange when possible, default: %d\n", m->re);
#ifdef IS2_MODEL
  fprintf(stderr, "  --nT:          set the number of temperatures, default: %d\n", m->ntp);
  fprintf(stderr, "  --T0:          set the minimal temperature, default: %g\n", m->tpmin);
  fprintf(stderr, "  --dT:          set the temperature increment, default: %g\n", m->tpdel);
  fprintf(stderr, "  --nsteps:      set the number of simulation steps, default: %d\n", m->nsteps);
  fprintf(stderr, "  --nequil:      set the number of equilibration steps, default: %d\n", m->nequil);
#endif
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc., default %d\n", m->verbose);
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
}



/* handle command line arguments */
__inline static void model_doargs(model_t *m, int argc, char **argv)
{
  int i, j, ch;
  char *p, *q;

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

      if ( strcmpfuzzy(p, "fnhis") == 0 ) {
        m->fnhis = q;
      } else if ( strcmpfuzzy(p, "fnhis2") == 0 ) {
        m->fnhis2 = q;
      } else if ( strcmpfuzzy(p, "de") == 0 ) {
        m->de = atof(q);
      } else if ( strcmpfuzzy(p, "dv") == 0 ) {
        m->dv = atof(q);
      } else if ( strcmpfuzzy(p, "wham") == 0 ) {
        m->wham_method = model_select(q, WHAM_NMETHODS, wham_methods);
      } else if ( strcmpfuzzy(p, "itmax") == 0 ) {
        m->itmax = atoi(q);
      } else if ( strcmpfuzzy(p, "tol") == 0 ) {
        m->tol = atof(q);
      } else if ( strcmpfuzzy(p, "nbases") == 0 ) {
        m->mdiis_nbases = atoi(q);
      } else if ( strcmpfuzzy(p, "mdamp") == 0 ) {
        m->mdiis_damp = atof(q);
      } else if ( strcmpfuzzy(p, "fndos") == 0 ) {
        m->fnlndos = q;
      } else if ( strcmpfuzzy(p, "fneav") == 0 ) {
        m->fneav = q;
      } else if ( strcmpfuzzy(p, "actmax") == 0 ) {
        m->actmax = atof(q);
      } else if ( strcmpfuzzy(p, "acmin") == 0 ) {
        m->acmin = atof(q);
      } else if ( strcmpfuzzy(p, "fnac") == 0 ) {
        m->fnac = q;
      } else if ( strcmpfuzzy(p, "re") == 0 ) {
        m->re = 1;
#ifdef IS2_MODEL
      } else if ( strcmpfuzzy(p, "nT") == 0 ) {
        m->ntp = atoi(q);
      } else if ( strcmpfuzzy(p, "T0") == 0 ) {
        m->tpmin = atof(q);
      } else if ( strcmpfuzzy(p, "dT") == 0 ) {
        m->tpdel = atof(q);
      } else if ( strcmpfuzzy(p, "nsteps") == 0 ) {
        m->nsteps = atoi(q);
      } else if ( strcmpfuzzy(p, "nequil") == 0 ) {
        m->nequil = atoi(q);
#endif
      } else if ( strcmpfuzzy(p, "help") == 0 ) {
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
      if ( strchr("ir", ch) != NULL ) {
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
        } else if ( ch == 'r' ) {
          m->radd = atof(q);
        }

        break; /* skip the rest of the characters in the option */
      } else if ( ch == 'H' ) {
        m->loadprev = 1;
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



#endif /* WHAMMODEL_H__ */
