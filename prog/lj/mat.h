#ifndef MAT_H__
#define MAT_H__



#include "vct.h"



/* a = b */
__inline static void mcopy(double a[D][D], double b[D][D])
{
  int i, j;

  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      a[i][j] = b[i][j];
    }
  }
}



/* a = b^T */
__inline static void mtrans(double a[D][D])
{
  int i, j;

  for ( i = 0; i < D; i++ ) {
    for ( j = i + 1; j < D; j++ ) {
      double x = a[i][j];
      a[i][j] = a[j][i];
      a[j][i] = x;
    }
  }
}



/* c = a^T b */
__inline static void mvtxv(double c[D][D], double a[D], double b[D])
{
  int i, j;

  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      c[i][j] = a[i] * b[j];
    }
  }
}



/* c = a b */
__inline static void mmxv(double c[D], double a[D][D], double b[D])
{
  int i, j, k;

  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      double s = 0;
      for ( k = 0; k < D; k++ )
        s += a[i][k] * b[k];
      c[i] = s;
    }
  }
}



/* c = a b */
__inline static void mmxm(double c[D][D], double a[D][D], double b[D][D])
{
  int i, j, k;

  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      double s = 0;
      for ( k = 0; k < D; k++ )
        s += a[i][k] * b[k][j];
      c[i][j] = s;
    }
  }
}



/* c = a^T b */
__inline static void mmtxm(double c[D][D], double a[D][D], double b[D][D])
{
  int i, j, k;

  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      double s = 0;
      for ( k = 0; k < D; k++ )
        s += a[k][i] * b[k][j];
      c[i][j] = s;
    }
  }
}



/* c = a b^T */
__inline static void mmxmt(double c[D][D], double a[D][D], double b[D][D])
{
  int i, j;

  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      c[i][j] = vdot(a[i], b[j]);
    }
  }
}



/* a += b * s */
__inline static void msinc(double a[D][D], double b[D][D], double s)
{
  int i, j;

  for ( i = 0; i < D; i++ ) {
    for ( j = 0; j < D; j++ ) {
      a[i][j] += b[i][j] * s;
    }
  }
}



/* compute the inverse matrix b = a^(-1), by Gaussian elimination */
__inline static int minv(double (*b)[D], double (*a)[D])
{
  int i, j, k, ip;
  double x;

  /* initialize b as the identity matrix */
  for ( i = 0; i < D; i++ )
    for ( j = 0; j < D; j++ )
      b[i][j] = (i == j);

  /* Gaussian elimination */
  for ( i = 0; i < D; i++ ) {
    /* choose the pivot as the largest element of column i */
    x = fabs( a[i][i] );
    for ( ip = i, k = ip + 1; k < D; k++ ) {
      if ( fabs( a[k][i] ) > x ) {
        ip = k;
        x = fabs( a[k][i] );
      }
    }

    /* swap the pivot (ip'th) row with the present row i */
    for ( k = i; k < D; k++ )
      x = a[i][k], a[i][k] = a[ip][k], a[ip][k] = x;
    for ( k = 0; k < D; k++ )
      x = b[i][k], b[i][k] = b[ip][k], b[ip][k] = x;

    /* normalize this row */
    x = a[i][i];
    if ( fabs(x) < DBL_MIN ) {
      fprintf(stderr, "Error: singular matrix\n");
      return -1;
    }
    for ( k = i; k < D; k++ ) a[i][k] /= x;
    for ( k = 0; k < D; k++ ) b[i][k] /= x;

    /* use the pivot row to zero the rest rows */
    for ( j = i + 1; j < D; j++ ) {
      x = a[j][i];
      for ( k = i; k < D; k++ )
        a[j][k] -= x * a[i][k];
      for ( k = 0; k < D; k++ )
        b[j][k] -= x * b[i][k];
    }
  }

  /* now that the matrix is upper triangular
   * make it diagonal */
  for ( i = D - 1; i >= 0; i-- ) {
    /* note a[i][i] should be 1 now */
    for ( j = 0; j < i; j++ ) {
      x = a[j][i];
      for ( k = 0; k < D; k++ )
        b[j][k] -= b[i][k] * x;
    }
  }
  return 0;
}




/* full pivot
 * return the pivot row r and column c, starting from (r0, c0)
 * cmap[r0] registers the actual column index */
__inline static double mpivotf_(double m[D][D], int r0,
    int cmap[], int *sgn)
{
  int i, j, r = r0, c = r0;
  double tmp, max, t;

  /* 1. find the pivot row and column */
  max = -1;
  for ( j = r0; j < D; j++ ) {
    for ( i = r0; i < D; i++ ) {
      if ( (tmp = fabs(m[i][j])) > max ) {
        r = i;
        c = j;
        max = tmp;
      }
    }
  }

  /* 2. put the pivot to the left-top corner */
  /* swap rows r and r0, which doesn't affect the solution */
  if ( r != r0 ) {
    vswap(m[r], m[r0]);
    if ( sgn ) *sgn *= -1;
  }

  if ( c != r0 ) { /* swap columns c and r0 */
    if ( sgn ) *sgn *= -1;
    for ( i = 0; i < D; i++ ) { /* must be from row 0 */
      t = m[i][c];
      m[i][c] = m[i][r0];
      m[i][r0] = t;
    }
    if ( cmap ) {
      i = cmap[c];
      cmap[c] = cmap[r0];
      cmap[r0] = i;
    }
  }

  return max;
}



/* return the determinant */
__inline static double mdet(double m[D][D])
{
  double y, det = 1.0, a[D][D];
  int i, j, k, sgn = 1;

  mcopy(a, m);
  for ( i = 0; i < D; i++ ) {
    /* find the pivot, the largest element in the matrix */
    if ( mpivotf_(a, i, NULL, &sgn) <= 0 )
      break;

    det *= a[i][i];

    /* normalize the row i */
    for ( y = 1.0 / a[i][i], k = i; k < D; k++ )
      a[i][k] *= y;

    /* use the pivot to simplify the matrix */
    for ( j = i + 1; j < D; j++ ) /* for rows j */
      for ( y = a[j][i], k = i; k < D; k++ ) /* for columns k >= i*/
        a[j][k] -= y * a[i][k];
  }

  return i < D ? 0 : det * sgn;
}



double msolvezero_lasty;
double msolvezero_lasttol;

/* Solve matrix equation a x = 0 by Gaussian elimination (full-pivot)
 * The matrix 'a' is destroyed, solutions are saved as *row* vectors in 'x'
 * return the number of solutions */
__inline static int msolvezero(double a[D][D], double (*x)[D], double reltol)
{
  double tol = 0, y;
  int i, j, k, cmap[D], sgn = 1;

  for ( i = 0; i < D; i++ ) {
    cmap[i] = i;
  }

  for ( i = 0; i < D; i++ ) {
    /* find the pivot, the largest element in the matrix */
    y = mpivotf_(a, i, cmap, &sgn);
    msolvezero_lasty = y;
    msolvezero_lasttol = tol;
    if ( y <= tol ) { /* we have D - i solutions */
      break;
    }
    if ( i == 0 ) {
      tol = y * reltol;
    }

    /* normalize the row i */
    y = 1.0 / a[i][i];
    for ( k = i; k < D; k++ ) {
      a[i][k] *= y;
    }

    /* use the pivot to simplify the matrix */
    for ( j = 0; j < D; j++ ) { /* for rows j */
      if ( j != i ) {
        for ( y = a[j][i], k = i; k < D; k++ ) { /* for columns k >= i*/
          a[j][k] -= y * a[i][k];
        }
      }
    }
  }

  /* solve the D - i solutions */
  for ( j = 0; j < D - i; j++ ) {
    vzero( x[j] );
    for ( k = 0; k < i; k++ ) {
      x[j][ cmap[k] ] = -a[k][i + j];
    }
    x[j][ cmap[i + j] ] = 1.0;
    vnormalize( x[j] );
  }

  return D - i;
}



/* given an eigenvalue, return the corresponding eigenvectors
 * Note: there might be multiple eigenvectors for the eigenvalue */
__inline static int meigvecs(double (*vecs)[D], double mat[D][D], double val)
{
  double m[D][D];
  int d;

  mcopy(m, mat); /* make a matrix */
  for ( d = 0; d < D; d++ ) {
    m[d][d] -= val;
  }
  return msolvezero(m, vecs, DBL_EPSILON * 10000.0);
}



/* sort `s' to descending order, order `u' and `v' correspondingly */
__inline static void msort2(double s[D], double (*u)[D], double (*v)[D])
{
  double t;
  int i, j, k;

  for ( i = 0; i < D; i ++ ) {
    for ( k = i, j = i + 1; j < D; j++ )
      if ( s[j] > s[k] )
        k = j;

    if ( k != i ) {
      t = s[i]; s[i] = s[k]; s[k] = t;
      if ( u ) {
        vswap( u[i], u[k] );
      }
      if ( v ) {
        vswap( v[i], v[k] );
      }
    }
  }
}



#if D == 3



/* compute eigenvalues of a 3x3 matrix
 * by solving a cubic equation */
__inline static double *meigval(double v[3], double a[3][3])
{
  double m, p, q, pr, pr3, a00, a11, a22;

  m = (a[0][0] + a[1][1] + a[2][2])/3;
  a00 = a[0][0] - m;
  a11 = a[1][1] - m;
  a22 = a[2][2] - m;
  q = ( a00 * (a11*a22 - a[1][2]*a[2][1])
      + a[0][1] * (a[1][2]*a[2][0] - a[1][0]*a22)
      + a[0][2] * (a[1][0]*a[2][1] - a11*a[2][0]) ) / 2.0;
  p = (a00*a00 + a11*a11 + a22*a22) / 6.0
    + (a[0][1]*a[1][0] + a[1][2]*a[2][1] + a[2][0]*a[0][2]) / 3.0;
  /* solve x^3 - 3 p x  - 2 q = 0 */
  pr = sqrt(p);
  pr3 = p * pr;
  if ( pr3 <= fabs(q) ) {
    if (q < 0.) { /* choose phi = pi/3 */
      v[1] = v[0] = m + pr;
      v[2] = m - 2.0 * pr;
    } else { /* phi = 0 */
      v[0] = m + 2.0 * pr;
      v[2] = v[1] = m - pr;
    }
  } else {
    double phi = acos(q/pr3)/3.0; /* 0 < phi < pi/3 */

    v[0] = m + 2.0 * pr * cos(phi);  /* largest */
    v[1] = m + 2.0 * pr * cos(phi - 2*M_PI/3); /* second largest */
    v[2] = m + 2.0 * pr * cos(phi + 2*M_PI/3); /* smallest */
  }
  return v;
}



/* given the matrix 'mat' and its eigenvalues 'v' return eigenvalues 'vecs'
 * ideally, eigenvalues should be sorted in magnitude-descending order
 * by default, vecs are transposed as a set of column vectors
 * set 'nt' != 0 to disable it: so vecs[0] is the first eigenvector  */
__inline static int meigsys(double v[3], double vecs[3][3], double mat[3][3], int nt)
{
  double vs[5][3] = {{0}}; /* for safety, vs needs 5 rows */
  int n = 0, nn, i = 0;

  meigval(v, mat);

  for ( nn = i = 0; i < 3; i++ ) {
    n = meigvecs(vs + nn, mat, v[nn]);
    if ( n == 0 ) {
      fprintf(stderr, "meigsys failed: try to increase msolvezero_reltol, i %d, nn %d, %g > %g\n",
          i, nn, msolvezero_lasty, msolvezero_lasttol);
      return -1;
    }
    if ( (nn += n) >= 3 ) break;
  }

  mcopy(vecs, vs);
  msort2(v, vecs, NULL);

  if ( !nt ) {
    mtrans(vecs);
  }
  return 0;
}



double msvd_reltol = 1e-6;

/* SVD decomposition of a matrix A = U S V^T */
__inline static void msvd(double a[3][3],
    double u[3][3], double s[3], double v[3][3])
{
  int i, rank;
  double ata[3][3], us[3][3];

  /* A^T A = V S^2 V^T, so (A^T A) V = V S^2 */

  /* 1. compute A^T A and its eigenvectors, which is V */
  mmtxm(ata, a, a);
  meigsys(s, v, ata, 1);

  /* 2. U^T = S^{-1} V^T A^T, and each row of U^T is an eigenvector
   * since eigenvectors are to be normalized, S^{-1} is unnecessary */
  if ( s[0] <= 0.0 ) {
    rank = 0;
    mcopy(u, v);
  } else {
    double tol = msvd_reltol;

    mmxmt(u, v, a);
    for ( i = 0; i < 3; i++ ) {
      vcopy(us[i], u[i]); /* save a copy of V^T A^T before normalizing it */
      s[i] = vnorm(u[i]);
      if ( s[i] > 0 ) {
        vsmul(u[i], 1 / s[i]);
      }
    }
    rank = 1;
    rank += (fabs( vdot(u[0], u[1]) ) < tol && s[1] > tol);
    rank += (fabs( vdot(u[0], u[2]) ) < tol
          && fabs( vdot(u[1], u[2]) ) < tol && s[2] > tol);
    if ( rank <= 2 ) {
      if ( rank == 1 ) {
        double z[3] = {0, 0, 0}, w, tmp;

        w = fabs( u[0][i = 0] );
        if ((tmp = fabs(u[0][1])) < w) w = tmp, i = 1;
        if ((tmp = fabs(u[0][2])) < w) i = 2;
        z[i] = 1.0f; /* select the smallest element in u[0] as z */
        vnormalize( vcross(u[1], z, u[0]) );
        s[1] = vdot(u[1], us[1]); /* S = U^T (V^T A^T)^T is more accurate than sqrt(A^T A) */
        if (s[1] < 0) { s[1] = -s[1]; vneg(u[1]); } /* make sure s[1] > 0 */
      }
      vnormalize( vcross(u[2], u[0], u[1]) );
      s[2] = vdot(u[2], us[2]);
      if (s[2] < 0) {
        s[2] = -s[2];
        vneg(u[2]);
      }
    }
    msort2(s, u, v);
  }
  mtrans(v);
  mtrans(u);
}



/* Fit x to y by rotation and translation of the `x'
 * If `refl', reflection can also be used.
 * The best-fit structure is saved to `xf', if not NULL */
__inline static double vrmsd(double (*x)[D], double (*xf)[D],
    double (*y)[D], const double *w, int n, int refl,
    double r[D][D], double *t)
{
  int i;
  double wi, wtot = 0, sq, dev = 0, dev0, detm;
  double xc[D], yc[D], xs[D], ys[D], xfi[D], sig[D], t_[D];
  double u[D][D], v[D][D], s[D][D] = {{0}}, xy[D][D], r_[D][D];

  if (r == NULL) r = r_;
  if (t == NULL) t = t_;

  /* 1. compute the centers */
  vzero(xc);
  vzero(yc);
  for ( wtot = 0., i = 0; i < n; i++ ) {
    wi = ( w != NULL ) ? w[i] : 1.0;
    vsinc(xc, x[i], wi);
    vsinc(yc, y[i], wi);
    wtot += wi;
  }
  vsmul(xc, 1.0/wtot);
  vsmul(yc, 1.0/wtot);

  /* 2. compute the asymmetric covariance matrix S = (x-xc) (y-yc)^T */
  for ( i = 0; i < n; i++ ) {
    wi = ( w != NULL ) ? w[i] : 1.0;

    vdiff(xs, x[i], xc); /* shift to the center avoid the translation */
    vdiff(ys, y[i], yc);
    mvtxv(xy, xs, ys);
    msinc(s, xy, wi);

    sq  = vsqr(xs);
    sq += vsqr(ys);
    dev += wi * sq; /* Tr(x^T x + y^T y) */
  }
  dev0 = dev;

  /* 3. SVD decompose S = u sig v^T */
  msvd(s, u, sig, v);

  /* 4. compute R = v u^T */
  mmxmt(r, v, u);
  detm = mdet(r);

  if ( detm < 0 && !refl ) { /* to avoid a reflection */
    mtrans(u);
    vneg(u[2]); /* flip the last eigenvector */
    mmxm(r, v, u);
    dev -= 2*(sig[0] + sig[1] - sig[2]);
    detm = mdet(r);
  } else {
    dev -= 2 * (sig[0] + sig[1] + sig[2]); /* -2 Tr(R x y^T) */
  }
  if ( dev < 0 ) {
    dev = 0;
  }
  mmxv(xs, r, xc); /* xs = R xc */
  vdiff(t, yc, xs); /* t = yc - R xc */

  /* 5. compute the rotated structure */
  if ( xf || dev < dev0 * 0.01 ) { /* if there's a large cancellation recompute the deviation */
    for ( dev = 0, i = 0; i < n; i++ ) {
      mmxv(xs, r, x[i]); /* xs = R x */
      vadd(xfi, xs, t); /* xfi = R x + t */
      sq = vdist2(y[i], xfi);
      if ( xf ) vcopy(xf[i], xfi);
      dev +=  (w ? w[i] * sq : sq); /* recompute the deviation */
    }
  }
  return sqrt(dev/wtot);
}



#endif /* D == 3 */



#endif /* MAT_H__ */

