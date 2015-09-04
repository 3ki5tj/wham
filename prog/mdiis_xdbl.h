#ifndef MDIIS_H__
#define MDIIS_H__



/* modified direct inversion of the iterative subspace (MDIIS) method
 * xdouble version */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "xdouble.h"
#include "lu_xdbl.h"



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
  int ibQ; /* index of the earliest base for the KTH scheme
              or that of the previous base for the HP scheme */
  xdouble (*getres)(void *, xdouble *, xdouble *); /* callback function */
  void *obj; /* object to pass to the callback function */
  xdouble **f;  /* basis */
  xdouble **res; /* residues */
  xdouble *mat; /* correlations of residues */
  xdouble *mat2; /* temporary matrix for LU decomposition */
  xdouble *coef; /* coefficients */
  xdouble *fbest;
  xdouble errmin;
  int verbose;
} mdiis_t;



/* open an mdiis object */
static mdiis_t *mdiis_open(int npt, int mnb,
    xdouble (*getres)(void *, xdouble *, xdouble *), void *obj,
    int verbose)
{
  mdiis_t *m;
  int mnb1;

  xnew(m, 1);
  m->npt = npt;
  m->mnb = mnb;
  m->nb = 0;
  m->ibQ = -1;
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
  for ( i = 0; i < nb; i++ ) {
    for ( j = 0; j < nb; j++ ) {
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
static void mdiis_gen(mdiis_t *m, xdouble *f,
    void (*normalize)(xdouble *, int, const void *),
    const void *obj, xdouble damp)
{
  int ib, il, npt = m->npt, nb = m->nb;

  for ( il = 0; il < npt; il++ )
    m->f[nb][il] = 0;
  for ( ib = 0; ib < nb; ib++ ) {
    xdouble coef = m->coef[ib];
    for ( il = 0; il < npt; il++ ) {
      xdouble x = m->f[ib][il] + damp * m->res[ib][il];
      m->f[nb][il] += coef * x;
    }
  }

  if ( normalize != NULL ) {
    normalize(m->f[nb], m->npt, obj);
  }

  /* f = m->f[nb] */
  cparr(f, m->f[nb], m->npt);
}



/* compute the dot product */
static xdouble mdiis_getdot(xdouble *a, xdouble *b, int n)
{
  int i;
  xdouble x = 0;

  for ( i = 0; i < n; i++ ) {
    x += a[i] * b[i];
  }
  return x / n;
}



/* build the residue correlation matrix */
static int mdiis_build(mdiis_t *m, xdouble *f, xdouble *res)
{
  int i, npt = m->npt;

  m->nb = 1;

  for ( i = 0; i < npt; i++ ) {
    m->f[0][i] = f[i];
    m->res[0][i] = res[i];
  }

  m->mat[0] = mdiis_getdot(m->res[0], m->res[0], npt);

  m->ibQ = -1;
  return 0;
}



/* try to add the new vector `f` and its residue `res`
 * into the base using the Kovalenko-Ten-no-Hirata scheme */
static int mdiis_update_kth(mdiis_t *m, xdouble *f, xdouble *res,
    xdouble err, xdouble threshold)
{
  int i, ibmin, ib, nb, mnb, npt = m->npt;
  xdouble dot, min;

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
    ib = (m->ibQ + 1) % m->mnb;
    m->ibQ = ib;
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
 * into the base using the Howard-Pettitt scheme */
static int mdiis_update_hp(mdiis_t *m, xdouble *f, xdouble *res,
    xdouble err)
{
  int i, ibmax, ib, nb, mnb, npt = m->npt;

  nb = m->nb;
  mnb = m->mnb;

  /* save this function if it achieves the minimal error so far */
  if ( err < m->errmin ) {
    cparr(m->fbest, m->f[nb], npt);
    m->errmin = err;
  }

  /* find the base with the largest residue */
  ibmax = 0;
  for ( i = 1; i < nb; i++ ) {
    /* the diagonal represents the error */
    if ( m->mat[i*mnb + i] > m->mat[ibmax*mnb + ibmax] ) {
      ibmax = i;
    }
  }

  /* if we are updating the same vector from the previous step
   * then we are stuck, so rebuild the basis
   * This condition cannot be true until we have a full basis */
  if ( ibmax == m->ibQ ) {
    /* rebuild the basis from f
     * if we rebuild the basis from f[ibmin]
     * it is more likely to enter a limit cycle */
    mdiis_build(m, f, res);
    return 0;
  }

  if ( nb < m->mnb ) {
    ib = nb;
    m->nb = ++nb;
  } else {
    ib = ibmax;
    /* Note: we only set ibQ if the basis is full */
    m->ibQ = ib;
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
 * into the base using the Howard-Pettitt scheme */
static int mdiis_update_hpl(mdiis_t *m, xdouble *f, xdouble *res,
    xdouble err)
{
  int i, ibmax, ib, nb, mnb, npt = m->npt;

  nb = m->nb;
  mnb = m->mnb;

  /* save this function if it achieves the minimal error so far */
  if ( err < m->errmin ) {
    cparr(m->fbest, m->f[nb], npt);
    m->errmin = err;
  }

  /* find the base with the largest (except the previous) residue */
  ibmax = ( m->ibQ == 0 && nb > 1 ) ? 1 : 0;
  for ( i = ibmax + 1; i < nb; i++ ) {
    if ( i == m->ibQ ) continue;
    /* the diagonal represents the error */
    if ( m->mat[i*mnb + i] > m->mat[ibmax*mnb + ibmax] ) {
      ibmax = i;
    }
  }

  if ( nb < m->mnb ) {
    ib = nb;
    m->nb = ++nb;
  } else {
    ib = ibmax;
    /* Note: we only set ibQ if the basis is full */
    m->ibQ = ib;
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
static int mdiis_update(mdiis_t *m, xdouble *f, xdouble *res,
    xdouble err)
{
  int i, ib, jb, nb, mnb, npt = m->npt;
  xdouble dot, max;

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

  if ( dot >= max ) {
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



__inline static void mdiis_summary(mdiis_t *m, int it)
{
  int j;
  fprintf(stderr, "it %d, %d bases:", it, m->nb);
  for ( j = 0; j < m->nb; j++ ) {
    fprintf(stderr, " %g(%g)", (double) m->coef[j],
        (double) SQRT(m->mat[j*m->mnb+j]));
  }
  fprintf(stderr, "\n");
}



__inline static int mdiis_log(mdiis_t *m,
    const xdouble *f, const xdouble *res,
    const xdouble *fref, int it)
{
  FILE *fp;
  char fn[FILENAME_MAX];
  int i, n = m->npt;

  sprintf(fn, "mdiis_iter%d.dat", it);
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write log %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d\n", it, n);
  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "%d %g %g", i, (double) f[i], (double) res[i]);
    if ( fref ) fprintf(fp, " %g", (double) (f[i] - fref[i]));
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}



enum {
  MDIIS_UPDATE_DEFAULT,
  MDIIS_UPDATE_KTH,
  MDIIS_UPDATE_HP,
  MDIIS_UPDATE_HPL,
  MDIIS_UPDATE_NMETHODS
};

const char *mdiis_update_methods[] = {
  "Default",
  "KTH",
  "HP",
  "HPL",
  "MDIIS_UPDATE_NMETHODS"
};



static xdouble iter_mdiis(xdouble *f, int npt,
    xdouble (*getres)(void *, xdouble *, xdouble *),
    void (*normalize)(xdouble *, int, const void *), void *obj,
    int nbases, xdouble damp,
    int update_method, xdouble threshold,
    int itmin, int itmax, xdouble tol,
    const xdouble *fref, int verbose)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib, success;
  xdouble err, errp, *res;
  clock_t t0, t1;

  t0 = clock();

  /* open an mdiis object */
  mdiis = mdiis_open(npt, nbases, getres, obj, verbose);
  /* use the space of the last array for the current residue */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial base set */
  mdiis->errmin = err = errp = mdiis->getres(obj, f, res);
  mdiis_build(mdiis, f, res);

  for ( it = 0; it < itmax; it++ ) {
    /* obtain a set of optimal coefficients of combination */
    mdiis_solve(mdiis);
    if ( verbose >= 2 ) {
      mdiis_summary(mdiis, it + 1);
    }

    /* generate a new f from the set of coefficients */
    mdiis_gen(mdiis, f, normalize, obj, damp);

    /* compute the residual vector `res` of `f` */
    err = mdiis->getres(obj, f, res);
    if ( fref != NULL ) {
      mdiis_log(mdiis, f, res, fref, it + 1);
    }

    /* add the new f into the basis */
    if ( update_method == MDIIS_UPDATE_KTH ) {
      ib = mdiis_update_kth(mdiis, f, res, err, threshold);
    } else if ( update_method == MDIIS_UPDATE_HP ) {
      ib = mdiis_update_hp(mdiis, f, res, err);
    } else if ( update_method == MDIIS_UPDATE_HPL ) {
      ib = mdiis_update_hpl(mdiis, f, res, err);
    } else {
      ib = mdiis_update(mdiis, f, res, err);
    }

    if ( verbose ) {
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d (%d)\n",
          it + 1, (double) errp, (double) err, ibp, ib, mdiis->nb);
    }
    if ( ib >= 0 ) {
      ibp = ib;
      errp = err;
    }

    if ( err < tol && it >= itmin ) {
      break;
    }
  }

  if ( err < tol ) {
    success = 1;
  } else { /* use the backup version */
    success = 0;
    cparr(f, mdiis->fbest, npt);
    err = mdiis->getres(obj, f, res);
  }

  mdiis_close(mdiis);

  t1 = clock();
  fprintf(stderr, "MDIIS finished in %d steps, error %g (%s), time %.4fs\n",
      it + 1, (double) err, (success ? "succeeded" : "failed"),
      1.0*(t1 - t0)/CLOCKS_PER_SEC);
  return err;
}



#endif /* MDIIS_H__ */

