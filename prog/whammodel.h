#ifndef WHAMMODEL_H__
#define WHAMMODEL_H__





/* This module mainly handles input parameters.
 * It is shared by multiple programs, and program-specific
 * parameters are embedded (which is not ideal). */





#include "util.h"





enum { SIMUL_MC, SIMUL_MD, SIMUL_NMETHODS };

const char *simul_methods[] = {"MC", "MD", "SIMUL_NMETHODS"};



typedef struct {
  char *prog;
  char *fnhis;
  char *fnhis2;
  int loadprev; /* load previous histogram */
  double radd; /* rate of adding trajectory frames into the histogram */
  int bootstrap;
  double de;
#ifndef IS2MODEL
  double dv;
#endif
#ifdef WHAM
  int wham_method;
#endif
#ifdef MBAR
  int mbar_method;
#endif
  int estimate;
  int itmin;
  int itmax;
  double tol;
  double damp;
  int mdiis_nbases; /* number of bases in MDIIS */
  int mdiis_update_method;
  double mdiis_threshold; /* controls when to clean up the basis in MDIIS */
  char *fnlndos;
  char *fneav;
  char *fnlndos2;
  char *fneav2;
  double actmax; /* cutoff time of the autocorrelation function */
  double acmin; /* minimum of the autocorrelation function */
  char *fnac;
  char *fnact;
  char *fninp;
  int verbose;
  int re;
#ifdef IS2_MODEL
  int L; /* override the length */
  int nT;
  double Tmin;
  double Tdel;
  int nequil;
  int nsteps;
#endif /* IS2_MODEL */
#ifdef LJ_MODEL
  int nn;
  double rho;
  double rcdef;
  double mcamp;
  double mddt;
  double thdt;
  double pdt;
  double emin;
  double emax;
  double vmin;
  double vmax;
  double Tmin;
  double Tdel;
  int nT;
  double Pmin;
  double Pdel;
  int nP;
  int nstadj;
  int nsteps;
  int nequil;
  int simul;
  double vamp;
  int nstvmov;
  int defsetup;
#endif /* LJ_MODEL */
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
  m->bootstrap = 0;
  m->de = 1.0;
  m->dv = 0.02;
#ifdef WHAM
  m->wham_method = WHAM_DIRECT;
#endif
#ifdef MBAR
  m->mbar_method = MBAR_DIRECT;
#endif
  m->estimate = 0;
  m->itmin = 0;
  m->itmax = 100000;
  m->tol = 1e-8;
  m->damp = 1.0;
  m->mdiis_nbases = 10;
#ifdef ENABLE_MDIIS
  m->mdiis_update_method = MDIIS_UPDATE_DEFAULT;
#endif
  m->mdiis_threshold = 10.0;
  m->fnlndos = NULL;
  m->fneav = NULL;
  m->fnlndos2 = NULL;
  m->fneav2 = NULL;
  m->actmax = 0;
  m->acmin = 0.05;
  m->fninp = NULL;
  m->fnac = NULL;
  m->fnact = NULL;
  m->verbose = 0;
  m->re = 0;
}



/* return the index of string from a predefined array */
__inline static int model_select(const char *s, int n, const char **arr)
{
  int i;

  for ( i = 0; i < n; i++ ) {
    if ( strncmpfuzzy(arr[i], s, strlen(arr[i])) == 0 ) {
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
  fprintf(stderr, "  -r:            set the rate of adding trajectory frames into the histogram, default %g\n", m->radd);
  fprintf(stderr, "  --bootstrap=:  set the number of tries for bootstrapping, default %d\n", m->bootstrap);
  fprintf(stderr, "  --de=:         set the energy bin size, default %g\n", m->de);
  fprintf(stderr, "  --dv=:         set the volume bin size, default %g\n", m->dv);
#ifdef WHAM
  fprintf(stderr, "  --wham=:       set the WHAM method, 'Direct', 'MDIIS', or 'ST', default: %s\n", wham_methods[m->wham_method]);
#endif
#ifdef MBAR
  fprintf(stderr, "  --mbar=:       set the MBAR method, 'Direct' or 'MDIIS', default: %s\n", mbar_methods[m->mbar_method]);
#endif
  fprintf(stderr, "  --est:         only estimate free energies, default %d\n", m->estimate);
  fprintf(stderr, "  --itmin=:      set the minimal number of iterations, default %d\n", m->itmin);
  fprintf(stderr, "  --itmax=:      set the maximal number of iterations, default %d\n", m->itmax);
  fprintf(stderr, "  --tol=:        set the tolerance of error, default %g\n", m->tol);
  fprintf(stderr, "  --damp=:       set the mixing factor in the direct method, default: %g\n", m->damp);
  fprintf(stderr, "  --nbases=:     set the number of bases in the MDIIS method, default: %d\n", m->mdiis_nbases);
  fprintf(stderr, "  --KTH:         use the Kovalenko-Ten-no-Hirata (queue-like) scheme to update the basis in the MDIIS method\n");
  fprintf(stderr, "  --HP:          use the Howard-Pettitt (dump the largest) scheme to update the basis in the MDIIS method\n");
  fprintf(stderr, "  --HPL:         use the modified Howard-Pettitt scheme to update the basis in the MDIIS method\n");
  fprintf(stderr, "  --mthreshold=: set the threshold to clean up the basis in the MDIIS method, default %g\n", m->mdiis_threshold);
  fprintf(stderr, "  --fndos=:      set the file for the density of states, default %s\n", m->fnlndos);
  fprintf(stderr, "  --fneav=:      set the file for the average energy, default %s\n", m->fneav);
  fprintf(stderr, "  --fndos2=:     set the file for the 2D density of states, default %s\n", m->fnlndos2);
  fprintf(stderr, "  --fneav2=:     set the file for the 2D average energy, default %s\n", m->fneav2);
  fprintf(stderr, "  --actmax=:     set the time cutoff of the autocorrelation function, default: %g\n", m->actmax);
  fprintf(stderr, "  --acmin=:      set thet minimal value of the autocorrelation function, default: %g\n", m->acmin);
  fprintf(stderr, "  --fnac=:       set the file for autocorrelation functions, default: %s\n", m->fnac);
  fprintf(stderr, "  --fnact=:      set the file for autocorrelation times, default: %s\n", m->fnact);
  fprintf(stderr, "  --re:          do replica exchange when possible, default: %d\n", m->re);
#ifdef IS2_MODEL
  fprintf(stderr, "  --L:           set the side length, default: %d\n", m->L);
  fprintf(stderr, "  --T0=:         set the minimal temperature, default: %g\n", m->Tmin);
  fprintf(stderr, "  --dT=:         set the temperature increment, default: %g\n", m->Tdel);
  fprintf(stderr, "  --nT=:         set the number of temperatures, default: %d\n", m->nT);
  fprintf(stderr, "  --nequil=:     set the number of equilibration steps, default: %d\n", m->nequil);
  fprintf(stderr, "  --nsteps=:     set the number of simulation steps, default: %d\n", m->nsteps);
#endif /* IS2_MODEL */
#ifdef LJ_MODEL
  fprintf(stderr, "  --nn=:         set the number of particles, default: %d\n", m->nn);
  fprintf(stderr, "  --rho=:        set the density, default: %g\n", m->rho);
  fprintf(stderr, "  --rcdef=:      set the preferred cutoff of the pair potential, default: %g\n", m->rcdef);
  fprintf(stderr, "  --mcamp:       set the Monte Carlo move size, default: %g\n", m->mcamp);
  fprintf(stderr, "  --mddt=:       set the MD time step, default: %g\n", m->mddt);
  fprintf(stderr, "  --thdt=:       set the thermostat time step, default: %g\n", m->thdt);
  fprintf(stderr, "  --pdt=:        set the barostat time step, default: %g\n", m->pdt);
  fprintf(stderr, "  --emin=:       set the minimal potential energy per particle, default: %g\n", m->emin);
  fprintf(stderr, "  --emax=:       set the maximal potential energy per particle, default: %g\n", m->emax);
  fprintf(stderr, "  --vmin=:       set the minimal volume per particle, default: %g\n", m->vmin);
  fprintf(stderr, "  --vmax=:       set the maximal volume per particle, default: %g\n", m->vmax);
  fprintf(stderr, "  --T0=:         set the minimal temperature, default: %g\n", m->Tmin);
  fprintf(stderr, "  --dT=:         set the temperature increment, default: %g\n", m->Tdel);
  fprintf(stderr, "  --nT=:         set the number of temperatures, default: %d\n", m->nT);
  fprintf(stderr, "  --P0=:         set the minimal pressure, default: %g\n", m->Pmin);
  fprintf(stderr, "  --dP=:         set the pressure increment, default: %g\n", m->Pdel);
  fprintf(stderr, "  --nP=:         set the number of pressure, default: %d\n", m->nP);
  fprintf(stderr, "  --nstadj=:     set the number of steps for adjustment, default: %d\n", m->nstadj);
  fprintf(stderr, "  --nequil=:     set the number of equilibration steps, default: %d\n", m->nequil);
  fprintf(stderr, "  --nsteps=:     set the number of simulation steps, default: %d\n", m->nsteps);
  fprintf(stderr, "  --mc:          do Monte Carlo simulations\n");
  fprintf(stderr, "  --md:          do molecular dynamics simulations\n");
  fprintf(stderr, "  --vamp=:       set the Monte Carlo log volume move size, default: %g\n", m->vamp);
  fprintf(stderr, "  --nstvmov=:    set the number of steps for volume moves, default: %d\n", m->nstvmov);
  fprintf(stderr, "  --defsetup:    use the default setup, default: %d\n", m->defsetup);
#endif /* LJ_MODEL */
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
      } else if ( strncmpfuzzy(p, "bootstrap", 4) == 0 ) {
        m->bootstrap = (q != NULL) ? atoi(q) : 1;
      } else if ( strcmpfuzzy(p, "de") == 0 ) {
        m->de = atof(q);
#ifndef IS2MODEL
      } else if ( strcmpfuzzy(p, "dv") == 0 ) {
        m->dv = atof(q);
#endif
#ifdef WHAM
      } else if ( strcmpfuzzy(p, "wham") == 0 ) {
        m->wham_method = model_select(q, WHAM_NMETHODS, wham_methods);
#endif
#ifdef MBAR
      } else if ( strcmpfuzzy(p, "mbar") == 0 ) {
        m->mbar_method = model_select(q, MBAR_NMETHODS, mbar_methods);
#endif
      } else if ( strcmpfuzzy(p, "est") == 0 ) {
        m->estimate = 1;
      } else if ( strcmpfuzzy(p, "itmin") == 0 ) {
        m->itmin = atoi(q);
      } else if ( strcmpfuzzy(p, "itmax") == 0 ) {
        m->itmax = atoi(q);
      } else if ( strcmpfuzzy(p, "tol") == 0 ) {
        m->tol = atof(q);
      } else if ( strcmpfuzzy(p, "damp") == 0
               || strcmpfuzzy(p, "mdamp") == 0 ) {
        m->damp = atof(q);
      } else if ( strcmpfuzzy(p, "nbases") == 0 ) {
        m->mdiis_nbases = atoi(q);
#ifdef ENABLE_MDIIS
      } else if ( strcmpfuzzy(p, "KTH") == 0 ) {
        m->mdiis_update_method = MDIIS_UPDATE_KTH;
      } else if ( strcmpfuzzy(p, "HP") == 0 ) {
        m->mdiis_update_method = MDIIS_UPDATE_HP;
      } else if ( strcmpfuzzy(p, "HPL") == 0 ) {
        m->mdiis_update_method = MDIIS_UPDATE_HPL;
#endif
      } else if ( strncmpfuzzy(p, "mthreshold", 4) == 0 ) {
        m->mdiis_threshold = atof(q);
      } else if ( strcmpfuzzy(p, "fndos") == 0 ) {
        m->fnlndos = q;
      } else if ( strcmpfuzzy(p, "fneav") == 0 ) {
        m->fneav = q;
      } else if ( strcmpfuzzy(p, "fndos2") == 0 ) {
        m->fnlndos2 = q;
      } else if ( strcmpfuzzy(p, "fneav2") == 0 ) {
        m->fneav2 = q;
      } else if ( strcmpfuzzy(p, "actmax") == 0 ) {
        m->actmax = atof(q);
      } else if ( strcmpfuzzy(p, "acmin") == 0 ) {
        m->acmin = atof(q);
      } else if ( strcmpfuzzy(p, "fnac") == 0 ) {
        m->fnac = q;
      } else if ( strcmpfuzzy(p, "fnact") == 0 ) {
        m->fnact = q;
      } else if ( strcmpfuzzy(p, "re") == 0 ) {
        m->re = 1;
#ifdef IS2_MODEL
      } else if ( strcmpfuzzy(p, "L") == 0 ) {
        m->L = atoi(q);
      } else if ( strcmpfuzzy(p, "nT") == 0 ) {
        m->nT = atoi(q);
      } else if ( strcmpfuzzy(p, "T0") == 0 ) {
        m->Tmin = atof(q);
      } else if ( strcmpfuzzy(p, "dT") == 0 ) {
        m->Tdel = atof(q);
      } else if ( strcmpfuzzy(p, "nequil") == 0 ) {
        m->nequil = atoi(q);
      } else if ( strcmpfuzzy(p, "nsteps") == 0 ) {
        m->nsteps = atoi(q);
#endif /* IS2_MODEL */
#ifdef LJ_MODEL
      } else if ( strcmpfuzzy(p, "nn") == 0 ) {
        m->nn = atoi(q);
      } else if ( strcmpfuzzy(p, "rho") == 0 ) {
        m->rho = atof(q);
      } else if ( strcmpfuzzy(p, "rcdef") == 0 ) {
        m->rcdef = atof(q);
      } else if ( strcmpfuzzy(p, "mcamp") == 0 ) {
        m->mcamp = atof(q);
      } else if ( strcmpfuzzy(p, "mddt") == 0 ) {
        m->mddt = atof(q);
      } else if ( strcmpfuzzy(p, "thdt") == 0 ) {
        m->thdt = atof(q);
      } else if ( strcmpfuzzy(p, "pdt") == 0 ) {
        m->pdt = atof(q);
      } else if ( strcmpfuzzy(p, "emin") == 0 ) {
        m->emin = atof(q);
      } else if ( strcmpfuzzy(p, "emax") == 0 ) {
        m->emax = atof(q);
      } else if ( strcmpfuzzy(p, "vmin") == 0 ) {
        m->vmin = atof(q);
      } else if ( strcmpfuzzy(p, "vmax") == 0 ) {
        m->vmax = atof(q);
      } else if ( strcmpfuzzy(p, "T0") == 0 ) {
        m->Tmin = atof(q);
      } else if ( strcmpfuzzy(p, "dT") == 0 ) {
        m->Tdel = atof(q);
      } else if ( strcmpfuzzy(p, "nT") == 0 ) {
        m->nT = atoi(q);
      } else if ( strcmpfuzzy(p, "P0") == 0 ) {
        m->Pmin = atof(q);
      } else if ( strcmpfuzzy(p, "dP") == 0 ) {
        m->Pdel = atof(q);
      } else if ( strcmpfuzzy(p, "nP") == 0 ) {
        m->nP = atoi(q);
      } else if ( strcmpfuzzy(p, "nstadj") == 0 ) {
        m->nstadj = atoi(q);
      } else if ( strcmpfuzzy(p, "nequil") == 0 ) {
        m->nequil = atoi(q);
      } else if ( strcmpfuzzy(p, "nsteps") == 0 ) {
        m->nsteps = atoi(q);
      } else if ( strcmpfuzzy(p, "mc") == 0 ) {
        m->simul = SIMUL_MC;
      } else if ( strcmpfuzzy(p, "md") == 0 ) {
        m->simul = SIMUL_MD;
      } else if ( strcmpfuzzy(p, "vamp") == 0 ) {
        m->vamp = atof(q);
      } else if ( strncmpfuzzy(p, "nstvmov", 7) == 0 ) {
        m->nstvmov = atoi(q);
      } else if ( strcmpfuzzy(p, "defsetup") == 0 ) {
        m->defsetup = 1;
#endif /* LJ_MODEL */
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
