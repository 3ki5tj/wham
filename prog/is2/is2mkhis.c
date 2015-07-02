/* Exact histograms for two-dimensional Ising model */
#ifndef IS2_LB
#define IS2_LB 5
#endif
#include "is2.h"
#define IS2_MODEL
#include "../whammodel.h"



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((size_t) (n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %ld\n", #x, (long) n); \
    exit(1); \
  }}
#endif



static void model_default_is2(model_t *m)
{
  model_default(m);
  m->de = 4;
  m->nT = 9;
  m->Tmin = 1.5;
  m->Tdel = 0.2;
  m->L = 1 << IS2_LB;
  m->fnlndos = "../../data/is2/islogdos32x32.txt";
  m->fnhis = "histref32x32.dat";
}



/* load density of states from file */
static int getlndos(double *arr, int n, const char *fn)
{
  FILE *fp;
  char s[1024], *p;
  int i;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open file %s\n", fn);
    return -1;
  }

  for ( i = 0; i <= n; i++ ) {
    fgets(s, sizeof s, fp);
    /* remove whatever that is after "`" */
    if ( ( p = strchr(s, '`')) != NULL )
      *p = '\0';
    sscanf(s, "%lf", &arr[i]);
  }

  return 0;
}



/* compute the histogram from the density of states */
static int mkhis(double **hs, double *lndos,
    double *beta, double *lnz, int nT, int n)
{
  int i, j;
  double e, x;

  for ( j = 0; j < nT; j++ ) {
    /* compute histogram j */
    for ( i = 0; i <= n; i++ ) {
      e = -2 * n + 4 * i;
      x = lndos[i] - beta[j] * e - lnz[j];
      hs[j][i] = exp(x);
    }
  }
  return 0;
}



/* save histograms to file */
static int savehis(double **hs, int nT, int n, const char *fn)
{
  int i, j;
  FILE *fp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write file %s\n", fn);
    return -1;
  }

  fprintf(fp, "# made by is2mkhis, n %d, nT %d\n", n, nT);

  for ( j = 0; j < nT; j++ ) {
    /* write histogram j */
    for ( i = 0; i <= n; i++ ) {
      if ( hs[j][i] < 1e-300 ) continue;
      fprintf(fp, "%d %g %d\n", -2 * n + i * 4, hs[j][i], j);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}



int main(int argc, char **argv)
{
  model_t m[1];
  double *beta, *lnz;
  double **hs, *lndos;
  int iT, n;

  model_default_is2(m);
  model_doargs(m, argc, argv);
  n = m->L * m->L;

  xnew(beta, m->nT);
  xnew(lnz, m->nT);
  xnew(hs, m->nT);
  for ( iT = 0; iT < m->nT; iT++ ) {
    beta[iT] = 1./(m->Tmin + m->Tdel * iT);
    lnz[iT] = is2_exact(m->L, m->L, beta[iT], NULL, NULL);
    xnew(hs[iT], n + 1);
  }
  xnew(lndos, n + 1);

  getlndos(lndos, n, m->fnlndos);
  mkhis(hs, lndos, beta, lnz, m->nT, n);
  savehis(hs, m->nT, n, m->fnhis);

  free(beta);
  free(lnz);
  for ( iT = 0; iT < m->nT; iT++ ) free(hs[iT]);
  free(hs);
  free(lndos);
  return 0;
}

