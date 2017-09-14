#ifndef VCT_H__
#define VCT_H__



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>



#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef D
#define D 3
#endif



typedef double vct[D];



__inline static double *vzero(double *x)
{
  int d;

  for ( d = 0; d < D; d++ )
    x[d] = 0;
  return x;
}



__inline static double *vneg(double *x)
{
  int d;

  for ( d = 0; d < D; d++ )
    x[d] = -x[d];
  return x;
}



__inline static double *vcopy(double *x, const double *y)
{
  int d;

  for ( d = 0; d < D; d++ )
    x[d] = y[d];
  return x;
}


__inline static void vswap(double *x, double *y)
{
  double z;
  int d;

  for ( d = 0; d < D; d++ ) {
    z = x[d];
    x[d] = y[d];
    y[d] = z;
  }
}



#define vinc(x, dx) vsinc(x, dx, 1)
#define vdec(x, dx) vsinc(x, dx, -1)

/* x += y * s */
__inline static double *vsinc(double *x, const double *y, double s)
{
  int d;

  for ( d = 0; d < D; d++ )
    x[d] += y[d] * s;
  return x;
}



#define vadd(c, a, b) vsadd(c, a, b, 1.0)
#define vdiff(c, a, b) vsadd(c, a, b, -1.0)

/* c = a + b * s */
__inline static double *vsadd(double *c, const double *a, const double *b, double s)
{
  int d;

  for ( d = 0; d < D; d++ )
    c[d] = a[d] + b[d] * s;
  return c;
}



/* c = -a - b */
__inline static double *vnadd(double *c, const double *a, const double *b)
{
  int d;

  for ( d = 0; d < D; d++ )
    c[d] = -a[d] - b[d];
  return c;
}



/* x *= s */
__inline static double *vsmul(double *x, double s)
{
  int d;

  for ( d = 0; d < D; d++ )
    x[d] *= s;
  return x;
}



/* y = x * s */
__inline static double *vsmul2(double *y, const double *x, double s)
{
  int d;

  for ( d = 0; d < D; d++ )
    y[d] = x[d] * s;
  return y;
}



#define vsqr(x) vdot(x, x)

__inline static double vdot(const double *x, const double *y)
{
  int d;
  double s = 0;

  for ( d = 0; d < D; d++ )
    s += x[d] * y[d];
  return s;
}



__inline static double *vwrap(double *x, double l)
{
  int d;

  for ( d = 0; d < D; d++ )
    x[d] = fmod(x[d] + 1000.*l, l);
  return x;
}



/* return the norm the vector */
__inline static double vnorm(double *a)
{
  return sqrt( vsqr(a) );
}



/* return the distance */
__inline static double vdistx(double *dx, const double *a, const double *b)
{
  return vnorm( vdiff(dx, a, b) );
}



/* return the square of the distance */
__inline static double vdist2(const double *a, const double *b)
{
  double dx[D];
  vdiff(dx, a, b);
  return vsqr( dx );
}



/* return the distance */
__inline static double vdist(const double *a, const double *b)
{
  double dx[D];
  return vnorm( vdiff(dx, a, b) );
}



/* normalize the vector */
__inline static double *vnormalize(double *v)
{
  double s = vsqr(v);
  s = ( s >= 0 ) ? 1./sqrt(s) : 1;
  return vsmul(v, s);
}



/* bond angle interaction */
__inline static double vang(const double *xi, const double *xj, const double *xk,
    double *gi, double *gj, double *gk)
{
  double xj_[D], xij[D], xkj[D], ri, rk, dot, ang;
  const double eps = DBL_EPSILON;

  if ( xj == NULL ) {
    /* use the origin as the default xj */
    vzero(xj_);
    xj = xj_;
  }

  ri = vdistx(xij, xi, xj);
  vsmul(xij, 1.0/ri);

  rk = vdistx(xkj, xk, xj);
  vsmul(xkj, 1.0/rk);

  dot = vdot(xij, xkj);
  if ( dot > 1.0 ) {
    dot = 1.0;
  } else if ( dot < -1.0 ) {
    dot = -1.0;
  }
  ang = acos( dot );

  if ( gi && gj && gk ) {
    double sn, gij, gkj;
    int d;
    sn = 1 - dot * dot;
    if ( sn < eps ) sn = eps;
    sn = -1.0 / sqrt(sn); /* -1.0/sin(phi) */
    for ( d = 0; d < D; d++ ) {
      gij = sn * (xkj[d] - xij[d]*dot) / ri;
      gkj = sn * (xij[d] - xkj[d]*dot) / rk;
      gi[d] = gij;
      gk[d] = gkj;
      gj[d] = -(gij + gkj);
    }
  }
  return ang;
}



#if D == 2


__inline static double vcross(double *x, double *y)
{
  return x[0]*y[1] - x[1]*y[0];
}



/* get a perpendicular vector to v */
__inline static double *vgetperp(double *p, const double *v)
{
  p[0] = -v[1];
  p[1] =  v[0];
  return p;
}


#elif D == 3



__inline static double *vcross(double *z, const double *x, const double *y)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
  return z;
}



/* get a perpendicular vector to v */
__inline static double *vgetperp(double *p, const double *v)
{
  int i, im = 0;
  double u[D], tmp, min;

  /* find the smallest component */
  min = fabs(v[im]);
  for ( i = 1; i < D; i++ ) {
    tmp = fabs(v[i]);
    if ( tmp < min ) {
      im = i;
      min = tmp;
    }
  }

  vzero(u);
  u[im] = 1;
  vnormalize(vcross(p, v, u));
  return p;
}


__inline static double vdih(const double *xi, const double *xj,
    const double *xk, const double *xl,
    double *gi, double *gj, double *gk, double *gl)
{
  double tol, phi, cosphi = 1;
  double nxkj, nxkj2, m2, n2;
  double xij[3], xkj[3], xkl[3], uvec[3], vvec[3], svec[3];
  double m[3], n[3]; /* the planar vector of xij x xkj,  and xkj x xkj */

  vdiff(xij, xi, xj);
  vdiff(xkj, xk, xj);
  vdiff(xkl, xk, xl);
  nxkj2 = vsqr(xkj);
  nxkj = sqrt(nxkj2);
  tol = nxkj2 * DBL_EPSILON;

  vcross(m, xij, xkj);
  m2 = vsqr(m);
  vcross(n, xkj, xkl);
  n2 = vsqr(n);
  if ( m2 > tol && n2 > tol ) {
    cosphi = vdot(m, n);
    cosphi /= sqrt(m2 * n2);
    if ( cosphi >= 1 ) {
      cosphi = 1;
    } else if ( cosphi < -1 ) {
      cosphi = -1;
    }
  }
  phi = acos(cosphi);
  if ( vdot(n, xij) < 0.0 ) phi = -phi;

  /* optionally calculate the gradient */
  if ( gi != NULL ) {
    if ( m2 > tol && n2 > tol ) {
      vsmul2(gi, m, nxkj/m2);
      vsmul2(gl, n, -nxkj/n2);
      vsmul2(uvec, gi, vdot(xij, xkj)/nxkj2);
      vsmul2(vvec, gl, vdot(xkl, xkj)/nxkj2);
      vdiff(svec, uvec, vvec);
      vdiff(gj, svec, gi);
      vnadd(gk, svec, gl);
    } else { /* clear the gradients */
      vzero(gi);
      vzero(gj);
      vzero(gk);
      vzero(gl);
    }
  }
  return phi;
}



#endif /* D == 3 */



#endif /* VCT_H__ */

