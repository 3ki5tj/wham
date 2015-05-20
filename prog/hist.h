#ifndef HIST_H__
#define HIST_H__



/* one-dimensional histograms */



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc((size_t) (n), sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %ld\n", #x, (long) n); \
    exit(1); \
  }}
#endif

#define HIST_VERBOSE    0x0001
#define HIST_ADDAHALF   0x0010
#define HIST_NOZEROES   0x0020
#define HIST_KEEPLEFT   0x0040
#define HIST_KEEPRIGHT  0x0080
#define HIST_KEEPLEFT2  0x0040
#define HIST_KEEPRIGHT2 0x0080
#define HIST_KEEPEDGE   (HIST_KEEPLEFT | HIST_KEEPRIGHT | HIST_KEEPLEFT2 | HIST_KEEPRIGHT2)
#define HIST_KEEPHIST   0x0100
#define HIST_OVERALL    0x0200
#define HIST_ADDITION   0x1000



typedef struct {
  int rows;
  int n;
  double xmin;
  double dx;
  double *arr;
} hist_t;



static void hist_clear(hist_t *hs)
{
  int i, n = hs->rows * hs->n;

  for ( i = 0; i < n; i++ )
    hs->arr[i] = 0;
}



static hist_t *hist_open(int rows, double xmin, double xmax, double dx)
{
  hist_t *hs;

  xnew(hs, 1);
  hs->rows = rows;
  hs->xmin = xmin;
  hs->dx   = dx;
  hs->n    = (int)((xmax - xmin)/dx + 0.99999999);
  xnew(hs->arr, hs->n * hs->rows);
  hist_clear(hs);
  return hs;
}



static void hist_close(hist_t *hs)
{
  free(hs->arr);
  free(hs);
}



/* compute sum, average and variance of the histogram h */
__inline static void hist_getsums(double *h, int n,
    double xmin, double dx, double *s)
{
  int i;
  double x, w;

  s[0] = s[1] = s[2] = 0.;
  for ( i = 0; i < n; i++ ) {
    x = xmin + (i + .5) * dx;
    w = h[i];
    s[0]  += w;
    s[1]  += w * x;
    s[2]  += w * x * x;
  }
  if ( s[0] > 0 ) {
    s[1] /= s[0];
    s[2] = s[2] / s[0] - s[1] * s[1];
  }
}



/* compute sum, average and variance of the rth histogram */
__inline static double hist_getave(const hist_t *hs, int r, double *s, double *var)
{
  double sums[3];

  hist_getsums(hs->arr + r * hs->n, hs->n, hs->xmin, hs->dx, sums);
  if ( s   != NULL ) *s   = sums[0];
  if ( var != NULL ) *var = sums[2];
  return sums[1];
}



/* write histograms to file */
__inline static int hist_save(const hist_t *hs, const char *fn, unsigned flags)
{
  const int version = 0;
  FILE *fp;
  int i, r, rp, rowp, imax, imin, rows = hs->rows, n = hs->n;
  const double *h = hs->arr, *p;
  double sm, (*sums)[3], fac, delta, xmin = hs->xmin, dx = hs->dx;
  double smtot[3] = {0}, *htot = NULL;

  if ((fp = fopen(fn, "w")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  /* get statistics */
  xnew(sums, rows);
  for ( r = 0; r < rows; r++ )
    hist_getsums(hs->arr + r * n, n, xmin, dx, sums[r]);

  /* compute the overall histogram */
  if (flags & HIST_OVERALL) {
    xnew(htot, n); /* the overall histogram */
    for (i = 0; i < n; i++) htot[i] = 0.;

    for (r = 0; r < rows; r++)
      for (i = 0; i < n; i++)
        htot[i] += h[r*n + i];

    hist_getsums(htot, n, xmin, dx, smtot);
    rowp = rows + 1;
  } else {
    rowp = rows;
  }

  /* print basic information */
  fprintf(fp, "# %d 0x%X | %d %d %g %g | ",
      version, flags, rows, n, xmin, dx);
  for (r = 0; r < rows; r++) /* number of visits */
    fprintf(fp, "%g ", sums[r][0]);
  fprintf(fp, "| ");
  for (r = 0; r < rows; r++) /* average, variance */
    fprintf(fp, "%g %g ", sums[r][1], sums[r][2]);
  fprintf(fp, "\n");

  delta = (flags & HIST_ADDAHALF) ? 0.5 : 0;

  for (r = 0; r < rowp; r++) {
    p = (r == rows) ? htot : (h+r*n);

    if (flags & HIST_KEEPRIGHT) {
      imax = n;
    } else { /* trim the right edge */
      for (i = n-1; i >= 0; i--)
        if (p[i] > 0)
          break;
      imax = i+1;
      if (imax == 0)
        continue;
    }

    if (flags & HIST_KEEPLEFT) {
      imin = 0;
    } else { /* trim the left edge */
      for (i = 0; i < imax; i++)
        if (p[i] > 0)
          break;
      imin = i;
    }

    sm = (r == rows) ? smtot[0] : sums[r][0];
    if (fabs(sm) < DBL_MIN) fac = 1.;
    else fac = 1.0/(sm*dx);

    for (i = imin; i < imax; i++) {
      if ((flags & HIST_NOZEROES) && p[i] < DBL_MIN)
        continue;
      fprintf(fp, "%g ", xmin+(i+delta)*dx);
      if (flags & HIST_KEEPHIST)
        fprintf(fp, "%.0f ", p[i]);
      rp = (r == rows) ? (-1) : r;
      fprintf(fp, "%20.14E %d\n", p[i]*fac, rp);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  if (flags & HIST_VERBOSE) {
    fprintf(stderr, "successfully wrote %s\n", fn);
    for (r = 0; r < rows; r++)
      fprintf(stderr, "%2d cnt: %20.4f av: %10.4f(%10.4f)\n",
          r, sums[r][0], sums[r][1], sums[r][2]);
  }
  free(sums);
  free(htot);
  return 0;
}



/* skip a | */
__inline static char *hist_skipabar(char *p)
{
  int next = -1;
  sscanf(p, " | %n", &next);
  return (next < 0) ? NULL : (p + next);
}



/* load a previous histogram
 * flags can have HIST_ADDITION and/or HIST_VERBOSE */
__inline static int hist_load(hist_t *hs, const char *fn, unsigned flags)
{
  FILE *fp;
  int verbose = (flags & HIST_VERBOSE);
  int add = (flags & HIST_ADDITION);
  int ver, next, hashist;
  int i, i1, r, r1, nlin = 0, rows = hs->rows, n = hs->n;
  unsigned fflags;
  double x, y, y2, fac, delta, xmin = hs->xmin, dx = hs->dx;
  double *sums = NULL;
  char s[65536] = "", *p;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }

  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  nlin++;
  if (6 != sscanf(s, " # %d 0x %X | %d%d%lf%lf | %n", &ver, &fflags, &r, &i, &y, &x, &next)
      || i < n || r != rows || fabs(x - dx) > 1e-5) {
    fprintf(stderr, "Error: bins = %d, %d, ng = %d, %d; dx = %g, %g\n",
        i, n, r, rows, x, dx);
    fclose(fp);
    return -1;
  }
  delta   = ((fflags & HIST_ADDAHALF) ? .5 : 0.);
  hashist =  (fflags & HIST_KEEPHIST);
  /* scan sums */
  xnew(sums, rows);
  for (p = s+next, r = 0; r < rows; r++) {
    if (1 != sscanf(p, "%lf%n", sums + r, &next)) {
      fprintf(stderr, "cannot read sums from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = hist_skipabar(p)) == NULL) {
    fprintf(stderr, "%s: no bar after the sums\n", fn);
    goto EXIT;
  }
  for (r = 0; r < rows; r++) {
    if (2 != sscanf(p, "%lf%lf%n", &y, &y2, &next)) {
      fprintf(stderr, "cannot read average/stddev from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }

  if (!add) { /* clear histogram */
    for (i = 0; i < rows * n; i++) hs->arr[i] = 0.;
  }

  /* loop over r = 0..rows-1 */
  for (r = 0; r < rows; r++) {
    fac = sums[r] * dx;
    while (fgets(s, sizeof s, fp)) {
      nlin++;
      for (p = s+strlen(s)-1; isspace((unsigned char)(*p)) && p >= s; p--)
        *p = '\0'; /* trim ending */
      if (s[0] == '#' || s[0] == '\0') break;
      if (hashist) {
        if (4 != sscanf(s, "%lf%lf%lf%d", &x, &y, &y2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      } else { /* simple */
        if (3 != sscanf(s, "%lf%lf%d", &x, &y2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      }
      if (r1 < 0) break; /* overall histogram */

      if (r1 < r) {
        fprintf(stderr, "wrong column index %d vs. %d on line %d, s=[%s]\n",
            r1, r, nlin, s);
        goto EXIT;
      } else if (r1 > r) {
        r = r1;
        fac = sums[r] * dx;
      }
      i1 = (int)((x - xmin)/dx - delta + .5);
      if (i1 < 0 || i1 >= n) {
        fprintf(stderr, "cannot find index for x = %g, delta = %g, i = %d/%d, on line %d\n",
            x, delta, i1, n, nlin);
        goto EXIT;
      }
      if ( !hashist ) y = y2*fac;
      if ( add ) y += hs->arr[r*n + i1];
      hs->arr[r*n + i1] = y;
    }
  }
  if (verbose)
    fprintf(stderr, "histogram loaded successfully from %s\n", fn);

  free(sums);
  fclose(fp);
  return 0;
EXIT:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  free(sums);
  /* clear histograms on error */
  for (i = 0; i < rows * n; i++) hs->arr[i] = 0.;
  return -1;
}



/* add x[r] of weight w into the rth row of the histogram, h, r = 0..rows-1
 * return the number of successful rows */
__inline static int histadd(const double *xarr, double w, double *h, int rows,
    int n, double xmin, double dx, unsigned flags)
{
  int r, ix, good = 0, verbose = flags & HIST_VERBOSE;
  double x;

  for (r = 0; r < rows; r++) {
    if ( (x = xarr[r]) < xmin ) {
      if (verbose)
        fprintf(stderr, "histadd underflows %d: %g < %g\n", r, x, xmin);
      continue;
    }
    if ( (ix = (int) ((x - xmin)/dx)) >= n ) {
      if (verbose)
        fprintf(stderr, "histadd overflows %d: %g > %g\n", r, x, xmin+dx*n);
      continue;
    }
    h[r*n + ix] += w;
    good++;
  }
  return good;
}



/* add x[r] of weight w into the rth row of the histogram */
__inline static int hist_add1(hist_t *hs, int r, double x, double w,
    unsigned flags)
{
  int ix, n = hs->n, verbose = flags & HIST_VERBOSE;
  double xmin = hs->xmin, dx = hs->dx;

  if ( r >= hs->rows || r < 0 ) {
    fprintf(stderr, "bad row index %d\n", r);
    return -1;
  }
  if ( x < xmin ) {
    if ( verbose )
      fprintf(stderr, "histadd underflows %d: %g < %g\n", r, x, xmin);
    return -1;
  }
  if ( (ix = (int) ((x - xmin)/dx)) >= n ) {
    if ( verbose )
      fprintf(stderr, "histadd overflows %d: %g > %g\n", r, x, xmin+dx*n);
    return -1;
  }
  hs->arr[r*n + ix] += w;
  return 0;
}



/* add x[r] of weight w into the rth histogram, r = 0..rows-1
 * return the number of successes */
__inline static int hist_add(hist_t *hs, const double *x, double w, unsigned flags)
{
  int r, good = 0;

  for ( r = 0; r < hs->rows; r++ )
    good += (hist_add1(hs, r, x[r], w, flags) == 0);
  return good;
}



/* fetch histogram size */
__inline static int hist_getinfo(const char *fn, int *row,
    double *xmin, double *xmax, double *xdel,
    int *version, unsigned *fflags)
{
  FILE *fp;
  char s[1024];
  int n;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  if (6 != sscanf(s, "# %d 0x %X | %d %d %lf %lf ",
        version, fflags, row, &n, xmin, xdel)) {
    fprintf(stderr, "%s: bad first line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  *xmax = *xmin + *xdel * n;
  fclose(fp);
  return 0;
}



/* initialize a histogram from file */
__inline static hist_t *hist_initf(const char *fn)
{
  int rows, version;
  unsigned fflags;
  double xmin, xmax, dx;
  hist_t *hs;

  if ( hist_getinfo(fn, &rows, &xmin, &xmax, &dx, &version, &fflags) != 0 ) {
    return NULL;
  }

  hs = hist_open(rows, xmin, xmax, dx);
  if ( hs == NULL ) {
    return NULL;
  }

  if ( hist_load(hs, fn, HIST_VERBOSE) != 0 ) {
    hist_close(hs);
    return NULL;
  }

  return hs;
}



#endif /* HIST_H__ */

