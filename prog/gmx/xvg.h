#ifndef XVG_H__
#define XVG_H__



/* loading data from the XVG file */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifdef MTRAND
#include "../mtrand.h"
#else
/* backup random number generator */
__inline static double rand01()
{
  static unsigned long s1 = 5318008, s2 = 11547776, s3 = 1040032;

#define TAUS(s,a,b,c,d) (((s &c) <<d) & 0xffffffffUL) ^ ((((s <<a) & 0xffffffffUL)^s) >>b)
  s1 = TAUS(s1, 13, 19, 4294967294UL, 12);
  s2 = TAUS(s2,  2, 25, 4294967288UL, 4);
  s3 = TAUS(s3,  3, 11, 4294967280UL, 17);
#undef TAUS
  return (s1^s2^s3)/4294967296.0;
}
#endif /* defined(MTRAND) */



const int xvg_blksz = 8192;

typedef struct {
  int ncap; /* capacity */
  int n;
  int m; /* number of y values */
  double *x;
  double **y;
  double dx;
} xvg_t;



__inline static xvg_t *xvg_open(int m)
{
  xvg_t *xvg;

  if ( (xvg = calloc(1, sizeof(*xvg))) == NULL ) {
    fprintf(stderr, "no memory for xvg\n");
    return NULL;
  }

  xvg->ncap = 0;
  xvg->n = 0;
  xvg->m = m;
  xvg->x = NULL;
  if ( (xvg->y = calloc(m, sizeof(*xvg->y))) == NULL ) {
    fprintf(stderr, "no memory for xvg->y\n");
    free(xvg);
    return NULL;
  }
  xvg->dx = 1;
  return xvg;
}



__inline static void xvg_close(xvg_t *xvg)
{
  int k;

  free(xvg->x);
  if ( xvg->y ) {
    for ( k = 0; k < xvg->m; k++ ) {
      free(xvg->y[k]);
    }
    free(xvg->y);
  }
  free(xvg);
}



/* detect the number of y rows */
__inline static int xvg_detectyrows(FILE *fp, double *dx)
{
  int m = 0;
  char buf[1024], *p, *q;
  double x1, x2;

  rewind(fp);
  /* get the first data line */
  while ( fgets(buf, sizeof(buf), fp) ) {
    if ( buf[0] != '#' && buf[0] != '@' ) {
      break;
    }
  }
  if ( feof(fp) ) {
    return 0;
  }

  /* compute the number of rows in buf */
  /* skip the leading spaces */
  for ( p = buf; *p && isspace(*p); p++ )
    ;

  sscanf(p, "%lf", &x1);

  while ( 1 ) {
    if ( *p == '\0' ) {
      break;
    }
    /* skip a token */
    for ( q = p; *q && !isspace(*q); q++ )
      ;
    m += 1;
    /* skip the space after the token */
    for ( p = q; *p && isspace(*p); p++ )
      ;
  }
  m -= 1;

  /* get the second data line */
  while ( fgets(buf, sizeof(buf), fp) ) {
    if ( buf[0] != '#' && buf[0] != '@' ) {
      break;
    }
  }
  if ( !feof(fp) ) {
    sscanf(buf, "%lf", &x2);
    *dx = x2 - x1;
  } else {
    *dx = 0;
  }

  return m;
}



/* build an xvg object from file */
#ifdef MTRAND
__inline static xvg_t *xvg_load(const char *fn, double r)
#else
__inline static xvg_t *xvg_load(const char *fn)
#endif
{
  xvg_t *xvg = NULL;
  int k, m, next;
  double xn, dx;
  char buf[1024], *p;
  FILE *fp;

  if ( fn == NULL || (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "xvg_load: cannot open %s\n", fn);
    return NULL;
  }

  /* detect the number of y rows */
  if ( (m = xvg_detectyrows(fp, &dx)) <= 0 ) {
    fprintf(stderr, "xvg_load: cannot determine the number of rows %d\n", m);
    return NULL;
  }

  if ( (xvg = xvg_open(m)) == NULL ) {
    return NULL;
  }

  xvg->dx = dx;

  rewind(fp);
  while ( fgets(buf, sizeof(buf), fp) ) {
    /* skip a comment line */
    if ( buf[0] == '#' || buf[0] == '@' ) {
      continue;
    }

#ifdef MTRAND
    if ( r < 1 && rand01() > r ) {
      continue;
    }
#endif

    /* reallocate memory */
    if ( xvg->n >= xvg->ncap ) {
      xvg->ncap += xvg_blksz;

      xvg->x = realloc(xvg->x, xvg->ncap * sizeof(xvg->x[0]));
      if ( xvg->x == NULL ) {
        fprintf(stderr, "no memory for xvg->x, %d\n", xvg->ncap);
        return NULL;
      }

      for ( k = 0; k < xvg->m; k++ ) {
        xvg->y[k] = realloc(xvg->y[k], xvg->ncap * sizeof(xvg->y[k][0]));
        if ( xvg->y[k] == NULL ) {
          fprintf(stderr, "no memory for xvg->y, %d\n", xvg->ncap);
          return NULL;
        }
      }
    }

    if ( 1 != sscanf(buf, "%lf%n", &xn, &next) ) {
      fprintf(stderr, "%s corrupted on line %d\n", fn, xvg->n);
      fclose(fp);
      xvg_close(xvg);
      return NULL;
    }

    /* if the time is not a multiple of the time step, drop it */
    if ( xvg->n >= 2
      && fmod(xn/xvg->dx + 1e-8, 1) > 1e-4 ) {
      //fprintf(stderr, "%s skip an intermediate step at %g, %g\n", fn, xn, fmod(xn/xvg->dx + 1e-8, 1));
      continue;
    }

    xvg->x[xvg->n] = xn;

    p = buf + next;
    for ( k = 0; k < xvg->m; k++ ) {
      if ( 1 != sscanf(p, "%lf%n", &xvg->y[k][xvg->n], &next) ) {
        fprintf(stderr, "%s corrupted on line %d\n", fn, xvg->n);
        fclose(fp);
        xvg_close(xvg);
        return NULL;
      }
      p += next;
    }
    xvg->n++;
  }

  fclose(fp);

  return xvg;
}



/* compute the minimal and maximal values */
__inline static void xvg_minmax(const xvg_t *xvg, double *ymin, double *ymax)
{
  int i, k, n = xvg->n, m = xvg->m;

  for ( k = 0; k < m; k++ ) { /* for m quantities */
    ymin[k] = 1e30;
    ymax[k] = -1e30;
    for ( i = 0; i < n; i++ ) {
      if ( xvg->y[k][i] < ymin[k] ) {
        ymin[k] = xvg->y[k][i];
      }
      if ( xvg->y[k][i] > ymax[k] ) {
        ymax[k] = xvg->y[k][i];
      }
    }
  }
}



/* compute the mean and variance */
__inline static void xvg_mean(const xvg_t *xvg, double *av)
{
  int i, k, n = xvg->n;
  double sy;

  for ( k = 0; k < xvg->m; k++ ) {
    sy = 0;
    for ( i = 0; i < n; i++ ) {
      sy += xvg->y[k][i];
    }
    av[k] = sy / n;
  }
}



/* construct a new trajectory from bootstrapping
 * `tau` is the autocorrelation time */
__inline static xvg_t *xvg_bootstrap(xvg_t* xvg0, double tau)
{
  xvg_t *xvg;
  int i, i0, k, n = xvg0->n, m = xvg0->m;
  double gam;

  xvg = xvg_open(m);
  xvg->n = n;
  xvg->ncap = n;
  xvg->x = calloc(n, sizeof(*xvg->x));
  xvg->dx = xvg0->dx;
  if ( xvg->x == NULL ) {
    fprintf(stderr, "no memory for xvg->x, %d\n", xvg->ncap);
    return NULL;
  }

  for ( k = 0; k < m; k++ ) {
    xvg->y[k] = calloc(n, sizeof(xvg->y[k][0]));
    if ( xvg->y[k] == NULL ) {
      fprintf(stderr, "no memory for xvg->y, %d\n", xvg->ncap);
      return NULL;
    }
  }

  /* bootstrapping */
  gam = (tau > 0 ? exp(-xvg->dx/tau) : 0);
  for ( i0 = i = 0; i < n; i++ ) {
    /* selective update the current frame */
    if ( rand01() > gam || i == 0 ) {
      i0 = (int) (rand01() * n);
    }
    xvg->x[i] = xvg0->x[i0];
    for ( k = 0; k < m; k++ ) {
      xvg->y[k][i] = xvg0->y[k][i0];
    }
  }
  return xvg;
}



/* compute the integral autocorrelation time `act[0 ... xvg->m-1]`
 * the maximal time for integration is given by `tcutoff`
 * the cutoff of the autocorrlation function is given by `min`
 * optionally save the autocorrelation function to `fn` */
__inline static int xvg_act(xvg_t *xvg, double *act,
    double tcutoff, double min, const char *fn)
{
  int i, j, k, done, n = xvg->n, m = xvg->m;
  double *av, *at;
  double dy1, dy2, syy, var = 0;
  FILE *fp = NULL;

  if ( (av = calloc(xvg->m, sizeof(*av))) == NULL ) {
    fprintf(stderr, "no memory for av\n");
    return -1;
  }

  if ( (at = calloc(xvg->m, sizeof(*at))) == NULL ) {
    fprintf(stderr, "no memory for at\n");
    return -1;
  }

  /* compute the mean */
  xvg_mean(xvg, av);

  for ( k = 0; k < m; k++ ) {
    act[k] = 0.5 * xvg->dx;
  }
  if ( fn != NULL ) {
    fp = fopen(fn, "w");
  }

  for ( j = 0; j < n - 1; j++ ) {
    if ( tcutoff > 0 && j * xvg->dx > tcutoff ) {
      break;
    }

    /* autocorrelation function at separation j */
    done = 1;
    for ( k = 0; k < m; k++ ) {
      syy = 0;
      for ( i = 0; i < n - j; i++ ) {
        dy1 = xvg->y[k][i] - av[k];
        dy2 = xvg->y[k][i+j] - av[k];
        syy += dy1 * dy2;
      }
      syy /= n - j;

      if ( j == 0 ) {
        var = syy;
        at[k] = 1.0;
        done = 0;
      } else {
        at[k] = syy / var;
        if ( at[k] >= min ) {
          act[k] += at[k] * xvg->dx;
          done = 0;
        }
      }
      //fprintf(stderr, "j %d, k %d, ac %g, tau %g\n", j, k, at[k], act[k]);
    }

    if ( fp != NULL ) {
      fprintf(fp, "%g", j * xvg->dx);
      for ( k = 0; k < m; k++ ) {
        fprintf(fp, " %g", at[k]);
      }
      fprintf(fp, "\n");
    }

    if ( done ) break;
  }

  if ( fp != NULL ) {
    fclose(fp);
  }

  free(av);
  free(at);
  return 0;
}



#endif /* XVG */
