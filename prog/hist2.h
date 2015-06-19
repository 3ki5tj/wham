#ifndef HIST2_H__
#define HIST2_H__



/* two-dimensional histograms */



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

#define HIST2_VERBOSE    0x0001
#define HIST2_ADDAHALF   0x0010
#define HIST2_NOZEROES   0x0020
#define HIST2_KEEPLEFT   0x0040
#define HIST2_KEEPRIGHT  0x0080
#define HIST2_KEEPLEFT2  0x0040
#define HIST2_KEEPRIGHT2 0x0080
#define HIST2_KEEPEDGE   (HIST2_KEEPLEFT | HIST2_KEEPRIGHT | HIST2_KEEPLEFT2 | HIST2_KEEPRIGHT2)
#define HIST2_KEEPHIST   0x0100
#define HIST2_OVERALL    0x0200
#define HIST2_INT        0x0400
#define HIST2_ADDITION   0x1000



typedef struct {
  int rows;
  int n, m;
  double xmin, ymin;
  double dx, dy;
  double *arr;
} hist2_t;



static void hist2_clear(hist2_t *hs2)
{
  int i, n = hs2->rows * hs2->n * hs2->m;
  for ( i = 0; i < n; i++ ) hs2->arr[i] = 0;
}



static hist2_t *hist2_open(int rows, double xmin, double xmax, double dx,
    double ymin, double ymax, double dy)
{
  hist2_t *hs2;

  xnew(hs2, 1);
  hs2->rows = rows;
  hs2->xmin = xmin;
  hs2->dx   = dx;
  hs2->n    = (int) ((xmax - xmin)/dx + 0.99999999);
  hs2->ymin = ymin;
  hs2->dy   = dy;
  hs2->m    = (int) ((ymax - ymin)/dy + 0.99999999);
  xnew(hs2->arr, hs2->n * hs2->m * hs2->rows);
  hist2_clear(hs2);
  return hs2;
}



static void hist2_close(hist2_t *hs2)
{
  free(hs2->arr);
  free(hs2);
}



__inline static void hist2_getsums(const double *h, int n,
    double xmin, double dx, int m, double ymin, double dy, double *s)
{
  double x, y, w;
  int i, j;

  s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0;
  for (i = 0; i < n; i++) {
    x = xmin + (i+.5)*dx;
    for (j = 0; j < m; j++) {
      y = ymin + (j+.5)*dy;
      w = h[i*m + j];
      s[0]  += w;
      s[1]  += w * x;
      s[2]  += w * y;
      s[3]  += w * x * x;
      s[4]  += w * x * y;
      s[5]  += w * y * y;
    }
  }
  if ( s[0] > 0 ) {
    s[1] /= s[0];
    s[2] /= s[0];
    s[3]  = s[3] / s[0] - s[1] * s[1];
    s[4]  = s[4] / s[0] - s[1] * s[2];
    s[5]  = s[5] / s[0] - s[2] * s[2];
  }
}



__inline static double hist2_getave(const hist2_t *hs, int r,
    double *s, double *avy, double *sxx, double *sxy, double *syy)
{
  int n = hs->n, m = hs->m;
  double sums[6];

  hist2_getsums(hs->arr + r * n * m,
      n, hs->xmin, hs->dx, m, hs->ymin, hs->dy, sums);
  if ( s      != NULL ) *s      = sums[0];
  if ( avy    != NULL ) *avy    = sums[2];
  if ( sxx    != NULL ) *sxx    = sums[3];
  if ( sxy    != NULL ) *sxy    = sums[4];
  if ( syy    != NULL ) *syy    = sums[5];
  return sums[1];
}



__inline static int hist2_save(const hist2_t *hs, const char *fn, unsigned flags)
{
  FILE *fp;
  int n = hs->n, m = hs->m, nm = n * m, rows = hs->rows;
  int i, j, r, imax, imin, jmax, jmin;
  double xmin = hs->xmin, dx = hs->dx, ymin = hs->ymin, dy = hs->dy;
  const double *p, *h = hs->arr;
  double (*sums)[6], fac, delta;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  nm = n * m;
  xnew(sums, rows);
  for ( r = 0; r < rows; r++ ) {
    hist2_getsums(hs->arr + r * nm, n, xmin, dx, m, ymin, dy, sums[r]);
  }
  /* print basic information */
  fprintf(fp, "# 1 0x%X | %d %d %g %g %d %g %g | ",
      flags, rows, n, xmin, dx, m, ymin, dy);
  for (r = 0; r < rows; r++) /* number of visits */
    fprintf(fp, "%20.14E ", sums[r][0]);
  fprintf(fp, " | ");
  for (r = 0; r < rows; r++) /* averages and standard deviations */
    fprintf(fp, "%g %g %g %g %g ", sums[r][1], sums[r][2],
        sums[r][3], sums[r][4], sums[r][5]);
  fprintf(fp, "| \n");

  delta = (flags & HIST2_ADDAHALF) ? 0.5 : 0;

  for (r = 0; r < rows; r++) { /* the rth data set */
    p = h + r * nm;

    if (flags & HIST2_KEEPRIGHT) {
      imax = n;
    } else { /* trim the right edge of i */
      for (i = n-1; i >= 0; i--) {
        for (j = 0; j < m; j++)
          if (p[i*m + j] > 0) break;
        if (j < m) break; /* found a nonzero entry */
      }
      imax = i + 1;
      if (imax == 0)
        continue;
    }

    if (flags & HIST2_KEEPLEFT) {
      imin = 0;
    } else { /* trim the left edge of i */
      for (i = 0; i < imax; i++) {
        for (j = 0; j < m; j++)
          if (p[i*m + j] > 0) break;
        if (j < m) break; /* found a nonzero entry */
      }
      imin = i;
    }

    if (flags & HIST2_KEEPRIGHT2) {
      jmax = m;
    } else { /* trim the right edge of j */
      for (j = m-1; j >= 0; j--) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmax = j + 1;
    }

    if (flags & HIST2_KEEPLEFT2) {
      jmin = 0;
    } else { /* trim the left edge of j */
      for (j = 0; j < jmax; j++) {
        for (i = imin; i < imax; i++)
          if (p[i*m + j] > 0) break;
        if (i < imax) break;
      }
      jmin = j;
    }

    if ( fabs(sums[r][0]) < 1e-6 ) fac = 1.;
    else fac = 1.0 / ( sums[r][0] * dx * dy );

    for (i = imin; i < imax; i++) {
      for (j = jmin; j < jmax; j++) {
        double x, y;
        if ((flags & HIST2_NOZEROES) && p[i*m + j] < 1e-16)
          continue;
        x = xmin + (i+delta)*dx;
        y = ymin + (j+delta)*dy;
        fprintf(fp, "%g %g ", x, y);
        if (flags & HIST2_KEEPHIST)
          fprintf(fp, "%20.14E ", p[i*m+j]);
        fprintf(fp, "%20.14E %d\n", p[i*m+j]*fac, r);
      }
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n#\n");
  }
  fclose(fp);
  if (flags & HIST2_VERBOSE) {
    fprintf(stderr, "successfully wrote %s\n", fn);
    for (r = 0; r < rows; r++) {
      double stdx = sqrt(sums[r][3]);
      double stdy = sqrt(sums[r][5]);
      fprintf(stderr, "%2d cnt: %20.4f xav: %10.4f(%10.4f) yav: %10.4f(%10.4f)\n",
          r, sums[r][0], sums[r][1], stdx, sums[r][2], stdy);
    }
  }
  free(sums);
  return 0;
}



/* skip a | */
__inline static char *hist2_skipabar(char *p)
{
  int next = -1;
  sscanf(p, " | %n", &next);
  return (next < 0) ? NULL : (p + next);
}



__inline static int hist2_load(hist2_t *hs, const char *fn, unsigned flags)
{
  FILE *fp;
  int verbose = (flags & HIST2_VERBOSE);
  int add = (flags & HIST2_ADDITION);
  int ver, next, hashist;
  int n = hs->n, m = hs->m, nm = n * m;
  int i, j, r, r1, nlin = 0, rows = hs->rows;
  unsigned fflags; /* flags from the file */
  double xmin = hs->xmin, dx = hs->dx, ymin = hs->ymin, dy = hs->dy;
  double x, y, xx, yy, xy, g, g2, fac, delta, *sums = NULL;
  double xmin1, dx1, ymin1, dy1;
  char s[40960] = "", *p;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }

  nm = n * m;
  /* check the first line */
  if (fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  nlin++;
  if (9 != sscanf(s, " # %d 0x %X | %d%d%lf%lf%d%lf%lf | %n", &ver, &fflags, &r,
        &i, &xmin1, &dx1, &j, &ymin1, &dy1, &next)
      || i < n || j < m || r != rows
      || fabs(dx1 - dx) > 1e-5 || fabs(dy1 - dy) > 1e-5 ) {
    fprintf(stderr, "%s error: bins %d, %d; %d, %d; ng %d, %d; dx %g, %g; dy %g, %g\n",
        fn, i, n, j, m, r, rows, dx1, dx, dy1, dy);
    fclose(fp);
    return -1;
  }
  delta   = ((fflags & HIST2_ADDAHALF) ? .5 : 0.);
  hashist =  (fflags & HIST2_KEEPHIST);
  /* scan sums */
  xnew(sums, rows);
  for (p = s+next, r = 0; r < rows; r++) {
    if (1 != sscanf(p, "%lf%n", sums + r, &next)) {
      fprintf(stderr, "cannot read sums from at %d/%d, s:\n%s\np:\n%s\n", r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = hist2_skipabar(p)) == NULL) goto EXIT;
  for (r = 0; r < rows; r++) {
    if (5 != sscanf(p, "%lf%lf%lf%lf%lf%n", &x, &y, &xx, &yy, &xy, &next)) {
      fprintf(stderr, "cannot read ave./cov. from at %d/%d, s:\n%s\np:\n%s\n",
          r, rows, s, p);
      goto EXIT;
    }
    p += next;
  }
  if ((p = hist2_skipabar(p)) == NULL) goto EXIT;

  if ( !add ) { /* clear histogram */
    for (i = 0; i < rows*nm; i++) hs->arr[i] = 0.;
  }

  /* loop over r = 0..rows-1 */
  for (r = 0; r < rows; r++) {
    fac = sums[r] * dx * dy;
    while ( fgets(s, sizeof s, fp) != NULL ) {
      nlin++;
      for ( p = s+strlen(s)-1; p >= s && isspace((unsigned char)(*p)); p-- )
        *p = '\0'; /* trim the ending */
      if (s[0] == '#') break;
      if (s[0] == '\0') continue;

      if (hashist) {
        if (5 != sscanf(s, "%lf%lf%lf%lf%d", &x, &y, &g, &g2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      } else {
        if (4 != sscanf(s, "%lf%lf%lf%d", &x, &y, &g2, &r1)) {
          fprintf(stderr, "error on line %d\n", nlin);
          goto EXIT;
        }
      }
      if ( r1 != r ) {
        fprintf(stderr, "wrong histogram index %d vs. %d on line %d\n",
          r1, r, nlin);
        goto EXIT;
      }
      i = (int)((x - xmin)/dx - delta + .5);
      if (i < 0 || i >= n) {
        fprintf(stderr, "cannot find index for x = %g\n", x);
        goto EXIT;
      }
      j = (int)((y - ymin)/dy - delta + .5);
      if (j < 0 || j >= m) {
        fprintf(stderr, "cannot find index for y = %g\n", y);
        return -1;
      }
      if ( !hashist ) g = g2 * fac;
      if ( flags & HIST2_INT ) {
        g = (double) (long) (g + 0.5);
      }
      if ( add ) g += hs->arr[r*nm + i*m + j];
      hs->arr[r*nm + i*m + j] = g;
    }
  }
  if (verbose) fprintf(stderr, "%s loaded successfully\n", fn);
  fclose(fp);
  return 0;
EXIT:
  fprintf(stderr, "error occurs at file %s, line %d, s:%s\n", fn, nlin, s);
  if (sums) free(sums);
  /* clear histograms on error */
  for (i = 0; i < rows * nm; i++) hs->arr[i] = 0.;
  fclose(fp);
  return -1;
}



/* add (x, y) into rth row of the histogram */
__inline static int hist2_add1(hist2_t *hs, int r,
    double x, double y, double w, unsigned flags)
{
  int ix, iy, n = hs->n, m = hs->m;
  int verbose = flags & HIST2_VERBOSE;
  double xmin = hs->xmin, dx = hs->dx, ymin = hs->ymin, dy = hs->dy;

  if ( x < xmin ) {
    if ( verbose )
      fprintf(stderr, "hist2add underflows, row %d: x %g < %g\n",
          r, x, xmin);
    return -1;
  }
  if ( y < ymin ) {
    if ( verbose )
      fprintf(stderr, "hist2add underflows, row %d: y %g < %g\n",
          r, y, ymin);
    return -1;
  }
  if ( (ix = (int)((x - xmin)/dx)) >= n ) {
    if (verbose)
      fprintf(stderr, "hist2add overflows, row %d: x %g > %g\n",
          r, x, xmin + dx*n);
    return -1;
  }
  if ( (iy = (int)((y - ymin)/dy)) >= m ) {
    if (verbose)
      fprintf(stderr, "hist2add overflows, row %d: y %g > %g\n",
          r, y, ymin + dy*m);
    return -1;
  }
  hs->arr[r*n*m + ix*m + iy] += w;
  return 0;
}



/* add (x[r], y[r]) of weight w into the rth histogram, r = 0..rows-1
 * return the number of successes */
__inline static int hist2_add(hist2_t *hs, const double *x, const double *y,
    int stride, double w, unsigned flags)
{
  int r, good = 0;

  for ( r = 0; r < hs->rows; r++ )
    good = (hist2_add1(hs, r, x[r * stride], y[r * stride], w, flags) == 0);
  return good;
}



/* fetch histogram information */
__inline static int hist2_getinfo(const char *fn, int *row,
    double *xmin, double *xmax, double *xdel,
    double *ymin, double *ymax, double *ydel,
    int *version, unsigned *fflags)
{
  FILE *fp;
  char s[1024];
  int n, m;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s: missing the first line\n", fn);
    fclose(fp);
    return -1;
  }
  if (9 != sscanf(s, "# %d 0x %X | %d %d %lf %lf %d %lf %lf ",
        version, fflags, row, &n, xmin, xdel, &m, ymin, ydel)) {
    fprintf(stderr, "%s: bad first line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  *xmax = *xmin + *xdel * n;
  *ymax = *ymin + *ydel * m;
  fclose(fp);
  return 0;
}



/* initialize a histogram from file */
__inline static hist2_t *hist2_initf(const char *fn, unsigned flags)
{
  int rows, version;
  unsigned fflags;
  double xmin, xmax, dx, ymin, ymax, dy;
  hist2_t *hs;

  if ( hist2_getinfo(fn, &rows, &xmin, &xmax, &dx,
        &ymin, &ymax, &dy, &version, &fflags) != 0 ) {
    return NULL;
  }

  hs = hist2_open(rows, xmin, xmax, dx, ymin, ymax, dy);
  if ( hs == NULL ) {
    return NULL;
  }

  if ( hist2_load(hs, fn, flags) != 0 ) {
    hist2_close(hs);
    return NULL;
  }

  return hs;
}



#endif /* HIST2_H__ */

