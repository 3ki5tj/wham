#ifndef XVG_H__
#define XVG_H__



/* loading data from the XVG file */



const int xvg_blksz = 1024;

typedef struct {
  int ncap; /* capacity */
  int n;
  double *x;
  double *y;
  double dx;
} xvg_t;



__inline static xvg_t *xvg_open(void)
{
  xvg_t *xvg;

  if ( (xvg = calloc(1, sizeof(*xvg))) == NULL ) {
    fprintf(stderr, "no memory for xvg\n");
    return NULL;
  }

  xvg->ncap = 0;
  xvg->n = 0;
  xvg->x = NULL;
  xvg->y = NULL;
  xvg->dx = 1;
  return xvg;
}



__inline static void xvg_close(xvg_t *xvg)
{
  free(xvg->x);
  free(xvg->y);
  free(xvg);
}



__inline static xvg_t *xvg_load(const char *fn)
{
  xvg_t *xvg = NULL;
  char buf[1024];
  FILE *fp;

  if ( fn == NULL || (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return xvg;
  }

  if ( (xvg = xvg_open()) == NULL ) {
    return xvg;
  }

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

      xvg->y = realloc(xvg->y, xvg->ncap * sizeof(*xvg->y));
      if ( xvg->y == NULL ) {
        fprintf(stderr, "no memory for y, %d\n", xvg->ncap);
        return NULL;
      }
    }

    sscanf(buf, "%lf %lf", &xvg->x[xvg->n], &xvg->y[xvg->n]);
    xvg->n++;
  }

  fclose(fp);

  if ( xvg->n >= 2 ) {
    xvg->dx = xvg->x[1] - xvg->x[0];
  }

  return xvg;
}



/* compute the minimal and maximal values */
__inline static double xvg_minmax(const xvg_t *xvg, double *ymax)
{
  int i, n = xvg->n;
  double ymin = 1e30;

  *ymax = -1e-30;
  for ( i = 0; i < n; i++ ) {
    if ( xvg->y[i] < ymin ) {
      ymin = xvg->y[i];
    }
    if ( xvg->y[i] > *ymax ) {
      *ymax = xvg->y[i];
    }
  }
  return ymin;
}



/* compute the mean and variance */
__inline static double xvg_mean(const xvg_t *xvg)
{
  int i, n = xvg->n;
  double sy;

  sy = 0;
  for ( i = 0; i < n; i++ ) {
    sy += xvg->y[i];
  }
  return sy / n;
}



/* compute the autocorrelation time
 * optionally save the autocorrelation function to `fn` */
__inline static double xvg_act(xvg_t *xvg,
    double tcutoff, double min, const char *fn)
{
  int i, k, n = xvg->n;
  double av, dy1, dy2, syy, var = 0, at, act;
  FILE *fp = NULL;

  av = xvg_mean(xvg);
  act = 0.5 * xvg->dx;
  if ( fn != NULL ) {
    fp = fopen(fn, "w");
  }

  for ( k = 0; k < n - 1; k++ ) {
    if ( tcutoff > 0 && k * xvg->dx > tcutoff ) {
      break;
    }

    /* autocorrelation function at separation k */
    syy = 0;
    for ( i = 0; i < n - k; i++ ) {
      dy1 = xvg->y[i] - av;
      dy2 = xvg->y[i+k] - av;
      syy += dy1 * dy2;
    }
    syy /= n - k;

    if ( k == 0 ) {
      var = syy;
    } else {
      at = syy / var;
      if ( at < min ) {
        break;
      }
      if ( fp != NULL ) {
        fprintf(fp, "%g %g\n", k * xvg->dx, at);
      }
      act += at * xvg->dx;
    }
  }

  if ( fp != NULL ) {
    fclose(fp);
  }
  return act;
}



#endif /* XVG */
