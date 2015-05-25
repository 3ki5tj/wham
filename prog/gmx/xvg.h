#ifndef XVG_H__
#define XVG_H__



/* loading data from the XVG file */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>



const int xvg_blksz = 1024;

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
__inline static int xvg_detectyrows(FILE *fp)
{
  int m = 0;
  char buf[1024], *p, *q;

  rewind(fp);
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
  return m;
}



/* build an xvg object from file */
__inline static xvg_t *xvg_load(const char *fn)
{
  xvg_t *xvg = NULL;
  int k, m, next;
  double xn;
  char buf[1024], *p;
  FILE *fp;

  if ( fn == NULL || (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "xvg_load: cannot open %s\n", fn);
    return xvg;
  }

  /* detect the number of y rows */
  if ( (m = xvg_detectyrows(fp)) <= 0 ) {
    return xvg;
  }

  if ( (xvg = xvg_open(m)) == NULL ) {
    return xvg;
  }

  rewind(fp);
  while ( fgets(buf, sizeof(buf), fp) ) {
    /* skip a comment line */
    if ( buf[0] == '#' || buf[0] == '@' ) {
      continue;
    }

    /* reallocate memory */
    if ( xvg->n >= xvg->ncap ) {
      xvg->ncap += xvg_blksz;

      xvg->x = realloc(xvg->x, xvg->ncap * sizeof(*xvg->x));
      if ( xvg->x == NULL ) {
        fprintf(stderr, "no memory for x, %d\n", xvg->ncap);
        return NULL;
      }

      for ( k = 0; k < xvg->m; k++ ) {
        xvg->y[k] = realloc(xvg->y[k], xvg->ncap * sizeof(xvg->y[k][0]));
        if ( xvg->y[k] == NULL ) {
          fprintf(stderr, "no memory for y, %d\n", xvg->ncap);
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
      //printf("%d %d %g\n", xvg->n, k, xvg->y[k][xvg->n]);
      p += next;
    }
    //printf("line %d\n", xvg->n); getchar();
    xvg->n++;
    if ( xvg->n == 2 ) {
      xvg->dx = xvg->x[1] - xvg->x[0];
    }
  }

  fclose(fp);

  return xvg;
}



/* compute the minimal and maximal values */
__inline static void xvg_minmax(const xvg_t *xvg, double *ymin, double *ymax)
{
  int i, k, n = xvg->n, m = xvg->m;

  for ( k = 0; k < m; k++ ) {
    ymin[k] = 1e30;
    ymax[k] = -1e30;
    for ( i = 0; i < n; i++ ) {
      //printf("k %d | i %d, y %g %g %g\n", k, i, xvg->y[k][i], ymin[k], ymax[k]);
      if ( xvg->y[k][i] < ymin[k] ) {
        ymin[k] = xvg->y[k][i];
      }
      if ( xvg->y[k][i] > ymax[k] ) {
        ymax[k] = xvg->y[k][i];
      }
    }
    //getchar();
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



/* compute the autocorrelation time
 * optionally save the autocorrelation function to `fn` */
__inline static int xvg_act(xvg_t *xvg, double *act,
    double tcutoff, double min, const char *fn)
{
  int i, j, k, n = xvg->n, m = xvg->m;
  double *av, *at;
  double dy1, dy2, syy, var = 0;
  FILE *fp = NULL;

  if ( (av = calloc(xvg->m, sizeof(*av))) == NULL ) {
    fprintf(stderr, "no memory for av\n");
    return -1;
  }

  if ( (at = calloc(xvg->m, sizeof(*act))) == NULL ) {
    fprintf(stderr, "no memory for act\n");
    return -1;
  }

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
      } else {
        at[k] = syy / var;
        if ( at[k] < min ) {
          break;
        }
        act[k] += at[k] * xvg->dx;
      }
    }

    if ( fp != NULL ) {
      fprintf(fp, "%g", j * xvg->dx);
      for ( k = 0; k < m; k++ ) {
        fprintf(fp, " %g", at[k]);
      }
      fprintf(fp, "\n");
    }
  }

  if ( fp != NULL ) {
    fclose(fp);
  }

  free(av);
  free(at);
  return 0;
}



#endif /* XVG */
