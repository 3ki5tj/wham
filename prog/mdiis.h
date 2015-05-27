#ifndef MDIIS_H__
#define MDIIS_H__



/* modified direct inversion of the iterative subspace (MDIIS) method */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "lu.h"



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %d\n", #x, (int) (n)); \
    exit(1); } }
#endif


#ifndef newarr
#define newarr(x, n) { int i_; xnew(x, n); \
  for ( i_ = 0; i_ < n; i_++ ) x[i_] = 0; }
#endif

/* copy array */
#ifndef cparr
#define cparr(x, y, n) { int i_; \
  for ( i_ = 0; i_ < n; i_++ ) x[i_] = y[i_]; }
#endif

#define delarr free

/* allocate a two-dimensional array */
#ifndef newarr2d
#define newarr2d(x, m, n) { int j_; \
  x = malloc(sizeof(x[0]) * m); \
  for ( j_ = 0; j_ < m; j_++ ) newarr( x[j_], n ); }
#endif

/* free a two-dimensional array */
#ifndef delarr2d
#define delarr2d(x, m) { int j_; \
  for ( j_ = 0; j_ < m; j_++ ) delarr( x[j_] ); \
  free(x); x = NULL; }
#endif





typedef struct {
  int npt;
  int mnb; /* maximal number of bases */
  int nb; /* number of functions in the basis */
  int ibQ; /* index of the earliest base */
  double (*getres)(void *, double *, double *); /* callback function */
  void *obj; /* object to pass to the callback function */
  double **f;  /* basis */
  double **res; /* residues */
  double *mat; /* correlations of residues */
  double *mat2; /* temporary matrix for LU decomposition */
  double *coef; /* coefficients */
  double *fbest;
  double errmin;
  int verbose;
} mdiis_t;



/* open an mdiis object */
static mdiis_t *mdiis_open(int npt, int mnb,
    double (*getres)(void *, double *, double *), void *obj,
    int verbose)
{
  mdiis_t *m;
  int mnb1;

  xnew(m, 1);
  m->npt = npt;
  m->mnb = mnb;
  m->nb = 0;
  m->ibQ = 0;
  m->getres = getres;
  m->obj = obj;
  mnb1 = mnb + 1;
  newarr2d(m->f,    mnb1, npt);
  newarr2d(m->res,  mnb1, npt);
  newarr(m->mat,    mnb * mnb);
  newarr(m->mat2,   mnb1 * mnb1);
  newarr(m->coef,   mnb1);
  newarr(m->fbest,  npt);
  m->errmin = DBL_MAX;
  m->verbose = verbose;
  return m;
}



/* close the mdiis object */
static void mdiis_close(mdiis_t *m)
{
  if ( m == NULL ) return;
  delarr2d(m->f,    m->mnb + 1);
  delarr2d(m->res,  m->mnb + 1);
  delarr(m->mat);
  delarr(m->mat2);
  delarr(m->coef);
  delarr(m->fbest);
  free(m);
}



/* solve the coefficients of combination */
static int mdiis_solve(mdiis_t *m)
{
  int nb = m->nb, nb1 = m->nb + 1, mnb = m->mnb, i, j;

  /* right-hand side of the equation */
  for ( i = 0; i < nb; i++ ) {
    m->coef[i] = 0;
  }
  m->coef[nb] = -1;

  /* copy the matrix, for the content is to be destroyed */
  for ( i = 0; i < nb1; i++ ) {
    for ( j = 0; j < nb1; j++ ) {
      m->mat2[i*nb1 + j] = m->mat[i*mnb + j];
    }
  }
  for ( i = 0; i < nb1; i++ ) {
    m->mat2[i*nb1 + nb] = -1;
    m->mat2[nb*nb1 + i] = -1;
  }
  m->mat2[nb*nb1 + nb] = 0;

  if ( lusolve(m->mat2, m->coef, nb1, 1e-20) != 0 ) {
    fprintf(stderr, "MDIIS lusolve failed\n");
    exit(1);
  }
  return 0;
}



/* construct the new f */
static void mdiis_gen(mdiis_t *m, double *f, double damp)
{
  int ib, il, npt = m->npt, nb = m->nb;

  for ( il = 0; il < npt; il++ )
    m->f[nb][il] = 0;
  for ( ib = 0; ib < nb; ib++ ) {
    double coef = m->coef[ib];
    for ( il = 0; il < npt; il++ ) {
      m->f[nb][il] += coef * (m->f[ib][il] + damp * m->res[ib][il]);
    }
  }

  /* f = m->f[nb] */
  cparr(f, m->f[nb], m->npt);
}



/* compute the dot product */
static double mdiis_getdot(double *a, double *b, int n)
{
  int i;
  double x = 0;

  for ( i = 0; i < n; i++ ) {
    x += a[i] * b[i];
  }
  return x / n;
}



/* build the residue correlation matrix */
static int mdiis_build(mdiis_t *m, double *f, double *res)
{
  int i, npt = m->npt;

  m->nb = 1;

  for ( i = 0; i < npt; i++ ) {
    m->f[0][i] = f[i];
    m->res[0][i] = res[i];
  }

  m->mat[0] = mdiis_getdot(m->res[0], m->res[0], npt);

  m->ibQ = 0;
  return 0;
}



/* try to add the new vector `f` and its residue `res`
 * into the base using the KTH scheme */
static int mdiis_update_kth(mdiis_t *m, double *f, double *res,
    double err, double threshold)
{
  int i, ibmin, ib, nb, mnb, npt = m->npt;
  double dot, min;

  nb = m->nb;
  mnb = m->mnb;

  /* save this function if it achieves the minimal error so far */
  if ( err < m->errmin ) {
    cparr(m->fbest, m->f[nb], npt);
    m->errmin = err;
  }

  /* choose the base with the smallest residue */
  ibmin = 0;
  for ( i = 1; i < nb; i++ ) {
    /* the diagonal represents the error */
    if ( m->mat[i*mnb + i] < m->mat[ibmin*mnb + ibmin] ) {
      ibmin = i;
    }
  }
  min = m->mat[ibmin*mnb + ibmin];
  dot = mdiis_getdot(res, res, npt);

  if ( dot > threshold * threshold * min ) {
    /* rebuild the basis */
    mdiis_build(m, m->f[ibmin], m->res[ibmin]);
    return 0;
  }

  if ( nb < m->mnb ) {
    ib = nb;
    m->nb = ++nb;
  } else {
    ib = m->ibQ; /* the earliest vector */
    m->ibQ = (ib + 1) % m->mnb;
  }

  /* replace base ib by f */
  for ( i = 0; i < npt; i++ ) {
    m->f[ib][i] = f[i];
    m->res[ib][i] = res[i];
  }

  /* update the residue correlation matrix
   * note: we do not need to update the last row & column */
  for ( i = 0; i < nb; i++ ) {
    m->mat[i*mnb + ib] = m->mat[ib*mnb + i]
      = mdiis_getdot(m->res[i], res, npt);
  }
  return ib;
}



/* try to add the new vector `f` and its residue `res`
 * into the base */
static int mdiis_update(mdiis_t *m, double *f, double *res,
    double err)
{
  int i, ib, jb, nb, mnb, npt = m->npt;
  double dot, max;

  nb = m->nb;
  mnb = m->mnb;

  /* save this function if it achieves the minimal error so far */
  if ( err < m->errmin ) {
    cparr(m->fbest, m->f[nb], npt);
    m->errmin = err;
  }

  /* choose the base with the largest residue */
  ib = 0;
  for ( i = 1; i < nb; i++ ) {
    /* the diagonal represents the error */
    if ( m->mat[i*mnb + i] > m->mat[ib*mnb + ib] ) {
      ib = i;
    }
  }
  max = m->mat[ib*mnb + ib];

  dot = mdiis_getdot(res, res, npt);

  if ( dot > max ) {
    if ( nb >= 2 ) {
      /* shrink the basis by removing the base with
       * the largest residue and try again */
      jb = nb - 1;
      if ( ib != jb ) {
        /* move the last base to position ib */
        for ( i = 0; i < npt; i++ ) {
          m->f[ib][i] = m->f[jb][i];
          m->res[ib][i] = m->res[jb][i];
        }

        for ( i = 0; i < nb - 1; i++ ) {
          if ( i == ib ) continue;
          m->mat[i*mnb + ib] = m->mat[i*mnb + jb];
          m->mat[ib*mnb + i] = m->mat[jb*mnb + i];
        }
        m->mat[ib*mnb + ib] = m->mat[jb*mnb + jb];
      }
      m->nb--;
      return -1;
    } else {
      /* rebuild the basis from f */
      mdiis_build(m, f, res);
      return 0;
    }
  }

  if ( nb < m->mnb ) {
    ib = nb;
    m->nb = ++nb;
  }

  /* replace base ib by f */
  for ( i = 0; i < npt; i++ ) {
    m->f[ib][i] = f[i];
    m->res[ib][i] = res[i];
  }

  /* update the residue correlation matrix
   * note: we do not need to update the last row & column */
  for ( i = 0; i < nb; i++ ) {
    m->mat[i*mnb + ib] = m->mat[ib*mnb + i]
      = mdiis_getdot(m->res[i], res, npt);
  }
  return ib;
}



static double iter_mdiis(double *f, int npt,
    double (*getres)(void *, double *, double *), void *obj,
    int nbases, double damp, int kth, double threshold,
    int itmax, double tol, int verbose)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib;
  double err, errp, *res;

  /* open an mdiis object */
  mdiis = mdiis_open(npt, nbases, getres, obj, verbose);
  /* use the space of the last array for the current residue */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial base set */
  mdiis->errmin = errp = mdiis->getres(obj, f, res);
  mdiis_build(mdiis, f, res);

  for ( it = 0; it < itmax && errp > tol; it++ ) {
    /* obtain a set of optimal coefficients of combination */
    mdiis_solve(mdiis);
    if ( verbose >= 2 ) {
      int j;
      fprintf(stderr, "it %d, %d bases:", it, mdiis->nb);
      for ( j = 0; j < mdiis->nb; j++ ) {
        fprintf(stderr, " %g(%g)", mdiis->coef[j], sqrt(mdiis->mat[j*mdiis->mnb+j]));
      }
      fprintf(stderr, "\n");
    }
    /* generate a new f from the set of coefficients */
    mdiis_gen(mdiis, f, damp);
    err = mdiis->getres(obj, f, res);
    /* add the new f into the basis */
    if ( kth ) {
      ib = mdiis_update_kth(mdiis, f, res, err, threshold);
    } else {
      ib = mdiis_update(mdiis, f, res, err);
    }

    if ( verbose ) {
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d (%d)\n",
          it, errp, err, ibp, ib, mdiis->nb);
    }
    if ( ib >= 0 ) {
      ibp = ib;
      errp = err;
    }
  }
  fprintf(stderr, "MDIIS finished in %d steps, err %g\n", it, errp);
  cparr(f, mdiis->fbest, npt);
  err = mdiis->getres(obj, f, res);
  mdiis_close(mdiis);
  return err;
}



#endif /* MDIIS_H__ */

