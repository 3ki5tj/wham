/* MBAR for a set of xvg files */
#include <time.h>
#include "../mtrand.h"
#define MTRAND
#include "xvg.h"
#define ENABLE_MDIIS
#include "mbar.h"
#define MBAR
#include "../whammodel.h"
#include "lsutil.h"





/* load trajectories */
static xvg_t **mkxvg(const char *fnls,
    int *nbeta, xdouble **beta, double radd)
{
  int i;
  char **fns;
  xvg_t **xvg = NULL;

  /* scramble the random number seed */
  mtscramble( time(NULL) + clock() );

  if ( (fns = getls(fnls, nbeta, beta)) == NULL ) {
    return NULL;
  }

  xnew(xvg, *nbeta);
  for ( i = 0; i < *nbeta; i++ ) {
    if ( (xvg[i] = xvg_load(fns[i], radd)) == NULL ) {
      exit(1);
    }
  }

  for ( i = 0; i < *nbeta; i++ ) {
    free(fns[i]);
  }
  free(fns);

  return xvg;
}



/* get the difference of lnz by BAR */
xdouble getdlnzbar(xvg_t *xvg0, xvg_t *xvg1, xdouble db,
    xdouble dlnz, xdouble tol, int itmax, int verbose)
{
  xdouble dlnn = LOG(1.0*xvg1->n/xvg1->n);
  xdouble e, lnden, lns0, lns1, res;
  xdouble dlnzp, resp;
  int it, i;

  /* exp(dlnz) =
   *      [ n0(E) + n1(E) ] exp(-db E)
   * Int ------------------------------- dE
   *       N0 + N1 exp(-db E -dlnz)
   * = N0 < exp(-db E) / [ N0 + N1 exp(-db E - dlnz) ] >_0
   * + N1 < exp(-db E) / [ N0 + N1 exp(-db E - dlnz) ] >_1
   *
   * exp(res)
   * = < 1 / [ exp(db E + dlnz) + N1/N0 ] >_0
   * + N1/N0 < 1 / [ exp(db E + dlnz) + N1/N0 ] >_1
   * */
  resp = 1e30;
  dlnzp = dlnz;
  for ( it = 0; it < itmax; it++ ) {
    lns0 = -1e30;
    for ( i = 0; i < xvg0->n; i++ ) {
      e = xvg0->y[0][i];
      lnden = mbar_lnadd(db * e + dlnz, dlnn);
      lns0 = mbar_lnadd(lns0, -lnden);
    }
    lns0 -= LOG(xvg0->n);

    lns1 = -1e30;
    for ( i = 0; i < xvg1->n; i++ ) {
      e = xvg1->y[0][i];
      lnden = mbar_lnadd(db * e + dlnz, dlnn);
      lns1 = mbar_lnadd(lns1, -lnden);
    }
    lns1 -= LOG(xvg1->n);

    res = mbar_lnadd(lns0, lns1 + dlnn);
    /* extrapolate to the point where res = 0
     * by the secant method */
    e = (resp * dlnz - res * dlnzp) / (resp - res);

    if ( verbose ) {
      fprintf(stderr, "it %d, res %g/%g, dlnz %.7f/%.7f/%.7f\n",
          it, (double) res, (double) resp, (double) dlnz, (double) dlnzp, (double) e);
    }
    if ( FABS(res) < tol ) break;

    resp = res;
    dlnzp = dlnz;
    if ( it == 0 ) {
      dlnz += res;
    } else {
      dlnz = e;
    }
  }

  return dlnz;
}



/* Gaussian partition: seek the position
 * where the energy distributions at two temperatures are the same */
static double getdlnzgp(double beta1, double eav1, double var1,
    double beta2, double eav2, double var2)
{
  double a1 = 1/var1, a2 = 1/var2;
  double A, B, C, D, num;

  A = a2 - a1;
  B = a2 * eav2 - a1 * eav1;
  /* the term log(a1/a2) is better dropped
   * for thermodynamic consistency */
  C = a2 * eav2 * eav2 - a1 * eav1 * eav1; /* + log(a1/a2); */
  D = sqrt(B*B - A*C);
  if ( eav1 > eav2 ) num = B + D;
  else num = B - D;
  return num/A * (beta1 - beta2);
}



/* Tilted Gaussian density of states */
static double getdlnztg(double beta1, double eav1, double var1,
    double beta2, double eav2, double var2)
{
  double de = eav2 - eav1, db = beta2 - beta1;
  return -(eav1 + eav2) / 2 * db + (1/var2 - 1/var1) * de * de / 12;
}



/* log variance */
static double getdlnzlnv(double beta1, double eav1, double var1,
    double beta2, double eav2, double var2)
{
  double de = eav2 - eav1, db = beta2 - beta1;
  return -(eav1 + eav2) / 2 * db + log(var2 / var1) * de * db / 12;
}



/* estimate free energies */
static int estimate(int nbeta, xvg_t **xvg, const xdouble *beta,
    xdouble tol, int itmax, int verbose)
{
  xdouble *eav, *var, *skw;
  xdouble *lnza, *lnzb, *lnzc, *lnzxpa, *lnzxpb, *lnzbar, *lnzgp, *lnztg, *lnzlnv;
  xdouble dlnza, dlnzb, dlnzc, dlnzxpa, dlnzxpb, dlnzbar, dlnzgp, dlnztg, dlnzlnv;
  xdouble emin, e, db;
  int i, j, n;

  xnew(eav, nbeta);
  xnew(var, nbeta);
  xnew(skw, nbeta);
  emin = 1e30;
  /* compute the variances */
  for ( j = 0; j < nbeta; j++ ) {
    eav[j] = var[j] = skw[j] = 0;
    n = xvg[j]->n;

    /* compute the mean */
    for ( i = 0; i < n; i++ ) {
      e = xvg[j]->y[0][i];
      eav[j] += e;
      if ( e < emin ) emin = e;
    }
    eav[j] /= n;

    /* compute the variance and skew */
    for ( i = 0; i < n; i++ ) {
      e = xvg[j]->y[0][i] - eav[j];
      var[j] += e * e;
      skw[j] += e * e * e;
    }
    if ( n > 1 ) {
      var[j] /= n - 1;
    }
    if ( n > 2 ) {
      skw[j] *= (xdouble) n / (n - 1) / (n - 2);
    }
    //printf("j %2d, beta %.6f, eav %.0f, var %8.0f, skew %8.0f\n", j, beta[j], eav[j], var[j], skw[j]);
    //if ( j > 0 ) printf("     db %.7f,             deav %8.0f, dvar %8.0f\n", beta[j]-beta[j-1], (eav[j]-eav[j-1])/(beta[j-1]-beta[j]), (var[j]-var[j-1])/(beta[j-1]-beta[j])); getchar();
  }

  emin -= 1;

  xnew(lnza, nbeta);
  xnew(lnzb, nbeta);
  xnew(lnzc, nbeta);
  xnew(lnzxpa, nbeta);
  xnew(lnzxpb, nbeta);
  xnew(lnzbar, nbeta);
  xnew(lnzgp, nbeta);
  xnew(lnztg, nbeta);
  xnew(lnzlnv, nbeta);

  lnza[0] = lnzb[0] = lnzc[0] = 0;
  lnzbar[0] = lnzgp[0] = lnztg[0] = lnzlnv[0] = 0;
  for ( j = 0; j < nbeta - 1; j++ ) {
    db = beta[j+1] - beta[j];
    dlnza = -(eav[j+1] + eav[j]) / 2 * db;
    dlnzb = dlnza - (var[j+1] - var[j]) / 12 * db * db;
    dlnzc = dlnza + (skw[j+1] + skw[j]) / 24 * db * db * db;
    lnza[j+1] = lnza[j] + dlnza;
    lnzb[j+1] = lnzb[j] + dlnzb;
    lnzc[j+1] = lnzc[j] + dlnzc;

    dlnzxpa = dlnzxpb = -1e30;

    /* from j to j+1 */
    n = xvg[j]->n;
    for ( i = 0; i < n; i++ ) {
      e = xvg[j]->y[0][i] - emin;
      dlnzxpa = mbar_lnadd(dlnzxpa, -db * e);
    }
    dlnzxpa += -db * emin - LOG(n);
    lnzxpa[j+1] = lnzxpa[j] + dlnzxpa;

    /* from j+1 to j */
    n = xvg[j+1]->n;
    for ( i = 0; i < n; i++ ) {
      e = xvg[j+1]->y[0][i] - emin;
      dlnzxpb = mbar_lnadd(dlnzxpb, db * e);
    }
    dlnzxpb += db * emin - LOG(n);
    lnzxpb[j+1] = lnzxpb[j] - dlnzxpb;

    dlnzbar = getdlnzbar(xvg[j], xvg[j+1], db, dlnzb,
        tol, itmax, verbose);
    lnzbar[j+1] = lnzbar[j] + dlnzbar;

    /* Gaussian partition, not so good */
    dlnzgp = getdlnzgp(beta[j], eav[j], var[j],
        beta[j+1], eav[j+1], var[j+1]);
    lnzgp[j+1] = lnzgp[j] + dlnzgp;

    /* tilted Gaussian, supposed to be good */
    dlnztg = getdlnztg(beta[j], eav[j], var[j],
        beta[j+1], eav[j+1], var[j+1]);
    lnztg[j+1] = lnztg[j] + dlnztg;

    /* similar to tilted the above, between dlnztg and dlnzb */
    dlnzlnv = getdlnzlnv(beta[j], eav[j], var[j],
        beta[j+1], eav[j+1], var[j+1]);
    lnzlnv[j+1] = lnzlnv[j] + dlnzlnv;
  }

  for ( j = 0; j < nbeta; j++ ) {
    printf("%3d %10.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f %15.7f %14.7f %14.7f %8d\n",
        j, (double) beta[j],
        (double) lnza[j], (double) lnzb[j], (double) lnzc[j],
        (double) lnzxpa[j], (double) lnzxpb[j], (double) lnzbar[j],
        (double) lnzgp[j], (double) lnztg[j], (double) lnzlnv[j],
        (double) eav[j], (double) SQRT(var[j]),
        (double) (BOLTZ*var[j]*beta[j]*beta[j]), xvg[j]->n);
  }

  free(eav);
  free(var);
  free(skw);
  free(lnza);
  free(lnzb);
  free(lnzc);
  free(lnzxpa);
  free(lnzxpb);
  free(lnzbar);
  free(lnzgp);
  free(lnztg);
  free(lnzlnv);
  return 0;
}



int main(int argc, char **argv)
{
  model_t m[1];
  xvg_t **xvg = NULL;
  int i, nbeta;
  unsigned flags = 0;
  xdouble *beta, *lnz;

  model_default(m);
  /* reduce the error tolerance for MBAR, if the precision is double */
  if ( strcmp(XDBLPRNF, "") == 0 ) {
    m->tol = 1e-7;
  }
  model_doargs(m, argc, argv);
  if ( m->rmcom ) flags |= MBAR_RMCOM;

  if ( m->fninp == NULL ) {
    model_help(m);
  }

  xvg = mkxvg(m->fninp, &nbeta, &beta, m->radd);
  if ( xvg == NULL ) {
    return -1;
  }

  if ( m->bootstrap ) {
    xvg_t *xvg0;
    double tau[10];

    /* scramble the random number seed */
    mtscramble( time(NULL) + clock() );
    /* bootstrapping */
    for ( i = 0; i < nbeta; i++ ) {
      xvg0 = xvg[i];
      xvg_act(xvg0, tau, m->actmax, m->acmin, m->fnac);
      xvg[i] = xvg_bootstrap( xvg0, tau[0] );
      xvg_close(xvg0);
    }
  }

  fprintf(stderr, "trajectories loaded, doing MBAR...\n");

  xnew(lnz, nbeta);
  for ( i = 0; i < nbeta; i++ ) {
    lnz[i] = 0;
  }

  if ( m->estimate ) {
    /* estimate */
    estimate(nbeta, xvg, beta,
        m->tol, m->itmax, m->verbose);
  } else {
    /* do MBAR */
    mbarx(nbeta, xvg, beta, lnz, flags,
        m->damp, m->mdiis_nbases,
        m->mdiis_update_method, m->mdiis_threshold,
        m->itmin, m->itmax, m->tol, m->verbose, m->mbar_method);

    if ( m->verbose ) {
      for ( i = 0; i < nbeta; i++ ) {
        printf("%3d %10.7f %14.7f %8d\n",
            i, (double) beta[i], (double) lnz[i], xvg[i]->n);
      }
    }
  }

  for ( i = 0; i < nbeta; i++ ) {
    free(xvg[i]);
  }
  free(xvg);
  free(beta);
  free(lnz);
  return 0;
}

