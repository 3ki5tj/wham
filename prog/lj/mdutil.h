#ifndef MDUTIL_H__
#define MDUTIL_H__



#include "mtrand.h"
#include "mat.h"



/* compute the center of mass */
static double getcom(double *xc, double (*x)[D], const double *m, int n)
{
  int i;
  double mtot = 0, wt;

  vzero(xc);
  for ( i = 0; i < n; i++ ) {
    wt = (m != NULL) ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);
  return mtot;
}

/* remove the center of mass motion */
static void rmcom(double (*x)[D], const double *m, int n)
{
  int i;
  double xc[D];

  getcom(xc, x, m, n);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], xc);
  }
}



#if D == 2



/* annihilate the total angular momentum by rotating the velocities */
__inline static void shiftangv(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  double u[D], xc[D], dx[D], newv;
  int i;

  /* compute the center of mass */
  getcom(xc, x, m, n);

  vzero(u);
  for ( i = 0; i < n; i++ ) {
    vdiff(dx, x[i], xc);
    /* u[0] and u[1] are the total angular momenta from
     * the current velocities, v[i], and the 90-degree
     * rotated velocities vi' = (-v[i][1], v[i][0])
     * and dx x vi' = dx[0]*v[i][0] - dx[1]*(-v[i][1]) */
    u[0] += m[i] * vcross(dx, v[i]);
    u[1] += m[i] * vdot(dx, v[i]);
  }
  vnormalize(u);

  /* if we now multiply the original velocity by u[1]
   * and multiply the rotated velocity by -u[0]
   * the total angular momentum would be zero
   * this is equivalent to a rotation with
   * cos(theta) = u[1] and sin(theta) = -u[0]
   * This choice of sign (instead of sin(theta) = -u[1], cos(theta) = u[0])
   * would ensure (v.x) is positive after the rotation */
  for ( i = 0; i < n; i++ ) {
    newv    =  u[1] * v[i][0] + u[0] * v[i][1];
    v[i][1] = -u[0] * v[i][0] + u[1] * v[i][1];
    v[i][0] = newv;
  }
}

/* annihilate the total angular momentum by rotation */
__inline static void shiftang(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i;
  double am, r2, wt, xc[D] = {0, 0}, xi[D];

  /* compute the center of mass */
  getcom(xc, x, m, n);

  am = r2 = 0.0;
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vdiff(xi, x[i], xc);
    am += wt * vcross(xi, v[i]);
    r2 += wt * vsqr(xi);
  }

  am = -am / r2;
  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    v[i][0] += -am * xi[1];
    v[i][1] +=  am * xi[0];
  }
}



#else



/* annihilate the total angular momentum */
__inline static void shiftangv(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i, j, k, round;
  double xc[D], xi[D], ang[D], am[D], z[D], th, rot[D][D];
  double wt, wdot, mat[D][D], vi[D], ama, aml;
  //double swd = 0;
  double tol = DBL_EPSILON * 1e3;

  /* compute the center of mass */
  getcom(xc, x, m, n);

  for ( round = 0; round < 100; round++ ) {
    /* compute the total angular momentum and
     * determine the axis of rotation, omega, such that
     *   Sum_i mi xi x vi is parallel to
     *   Sum_i mi xi x (omega x vi)
     *   = Sum_i mi [(vi.xi) - (vi::xi)] omega */
    vzero(am);
    mzero(mat);
    ama = 0;
    for ( i = 0; i < n; i++ ) {
      wt = ( m != NULL ) ? m[i] : 1.0;
      vdiff(xi, x[i], xc);
      vcross(ang, xi, v[i]);
      vsinc(am, ang, wt);
      /* this sum is used to estimate the magnitude of the
       * the angular momentum to establish an error threshold */
      ama += wt * vnorm(ang);
      wdot = wt * vdot(v[i], xi);
      for ( j = 0; j < D; j++ ) {
        mat[j][j] += wdot;
        for ( k = 0; k < D; k++ ) {
          mat[j][k] -= wt * v[i][j] * xi[k];
        }
      }
      //swd += wdot;
    }
    aml = vnorm(am);
    if ( aml < ama * tol ) break;
    msolve(mat, am);
    /* now the direction of `am` reprensent the axis of rotation,
     * omega, around which an 90-degree rotation (omega x vi) with
     * a scaling of |am| will produce the same angular momentum
     * as the initial velocity, vi.
     * which means the ratio of the old and new velocity field
     * should be 1 to -|am| */
    th = atan2(-vnorm(am), 1);
    //printf("round %d, th %g, am %g, ama %g, |L| %g, x.v %g\n", round, th, vnorm(am), ama, aml, swd);
    vnormalize(vcopy(z, am));

    mrota(rot, z, th);

    for ( i = 0; i < n; i++ ) {
      mmxv(vi, rot, v[i]);
      vcopy(v[i], vi);
    }
  }

  /* additional 180-degree to maximize mi(ri.vi) */
  {
    double val[D], tr;
    mzero(mat);
    for ( i = 0; i < n; i++ ) {
      wt = ( m != NULL ) ? m[i] : 1.0;
      vdiff(xi, x[i], xc);
      for ( j = 0; j < D; j++ ) {
        for ( k = 0; k < D; k++ ) {
          mat[j][k] += wt * xi[j] * v[i][k];
        }
      }
    }
    for ( tr = 0, j = 0; j < D; j++ ) tr += mat[j][j];
    meigsys(val, rot, mat, 1);
    if ( val[0] > tr ) {
      /* rotate around rot[0] for 180 degrees */
      for ( i = 0; i < n; i++ ) {
        double dot = vdot(v[i], rot[0]);
        vsinc(v[i], rot[0], -2*dot);
        vneg(v[i]);
      }
    }
  }
}



/* annihilate the total angular momentum
 * solve
 *   /  m (y^2 + z^2)   -m x y          -m x y        \
 *   |  -m x y          m (x^2 + z^2)   -m y z        |  c  =  L
 *   \  -m x z          -m y z          m (x^2 + y^2) /
 * use a velocity field
 *    v' = v - c x r
 *   */
__inline static void shiftang(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i;
  double xc[D], xi[D], ang[D], am[D];
  double dv[D], mat[D][D];
  double xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;
  double wt;

  /* compute the center of mass */
  getcom(xc, x, m, n);

  vzero(am);
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vdiff(xi, x[i], xc);
    vcross(ang, xi, v[i]);
    vsinc(am, ang, wt);

    xx += wt * xi[0] * xi[0];
    yy += wt * xi[1] * xi[1];
    zz += wt * xi[2] * xi[2];
    xy += wt * xi[0] * xi[1];
    yz += wt * xi[1] * xi[2];
    zx += wt * xi[2] * xi[0];
  }
  mat[0][0] = yy + zz;
  mat[1][1] = xx + zz;
  mat[2][2] = xx + yy;
  mat[0][1] = mat[1][0] = -xy;
  mat[1][2] = mat[2][1] = -yz;
  mat[0][2] = mat[2][0] = -zx;

  /* the axis of rotation is the solution of -M^(-1) * L */
  msolve(mat, am);
  vneg(am);

  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    vcross(dv, am, xi);
    vinc(v[i], dv);
  }
}



#endif /* D == 3 */



/* compute the kinetic energy */
__inline static double md_ekin(double (*v)[D], const double *m, int n)
{
  int i;
  double ek = 0, wt;

  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1;
    ek += wt * vsqr( v[i] );
  }
  return ek * 0.5;
}



/* exact velocity rescaling thermostat */
__inline static double md_vrescale(double (*v)[D],
    const double *m, int n, int dof, double tp, double dt)
{
  int i;
  double ek1, ek2, s, c, r, r2;

  c = exp(-dt);
  ek1 = md_ekin(v, m, n);
  r = randgaus();
  r2 = randchisqr(dof - 1);
  ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * tp / 2);
  if ( ek2 < 0 ) {
    ek2 = 0;
  }
  s = sqrt(ek2 / ek1);
  for (i = 0; i < n; i++) {
    vsmul(v[i], s);
  }
  return ek2;
}



/* Nose-Hoover chain thermostat */
__inline static double md_nhchain(double (*v)[D],
    const double *m, int n, int dof, double tp, double dt,
    int nnhc, double *zeta, const double *zmass)
{
  int i, j;
  double s, GQ, mvv;

  mvv = md_ekin(v, m, n) * 2;

  for ( j = nnhc - 1; j >= 0; j-- ) {
    s = ( j == nnhc - 1 ) ? 1 : exp(-zeta[j+1]*dt*0.25);
    GQ = ( j == 0 ) ? (mvv - dof * tp) : (zmass[j-1] * zeta[j-1] * zeta[j-1] - tp);
    zeta[j] = (zeta[j] * s + GQ /zmass[j] * dt*0.5) * s;
  }

  s = exp( -zeta[0] * dt );
  for ( i = 0; i < n; i++ ) {
    vsmul(v[i], s);
  }
  mvv *= s * s;

  for ( j = 0; j < nnhc; j++ ) {
    s = ( j == nnhc - 1 ) ? 1 : exp(-zeta[j+1]*dt*0.25);
    GQ = ( j == 0 ) ? (mvv - dof * tp) : (zmass[j-1] * zeta[j-1] * zeta[j-1] - tp);
    zeta[j] = (zeta[j] * s + GQ /zmass[j] * dt*0.5) * s;
  }

  return mvv * 0.5;
}




__inline static double md_langevin(double (*v)[D],
    const double *m, int n, double tp, double dt)
{
  int i, k;
  double s, v0, vi, wt, ek = 0;

  s = exp(-dt);
  v0 = sqrt( tp * (1 - s * s) );
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1;
    vi = v0 / sqrt( wt );
    for ( k = 0; k < D; k++ ) {
      v[i][k] = v[i][k] * s + vi * randgaus();
    }
    ek += 0.5 * wt * vsqr(v[i]);
  }
  return ek;
}



/* randomly swap the velocities of k pairs of particles */
__inline static double md_vscramble(double (*v)[D],
    const double *m, int n, int k)
{
  int i, j, l, d;
  double vi, vj, smi, smj;

  for ( l = 0; l < k; l++ ) {
    i = (int) (rand01() * n);
    j = (i + 1 + (int) (rand01() * (n - 1))) % n;

    if ( m != NULL ) {
      smi = sqrt( m[i] );
      smj = sqrt( m[j] );
    } else {
      smi = smj = 1.0;
    }

    for ( d = 0; d < D; d++ ) {
      vi = smi * v[i][d];
      vj = smj * v[j][d];
      v[i][d] = vj / smi;
      v[j][d] = vi / smj;
    }
  }
  return md_ekin(v, m, n);
}



/* bond energy k (r - r0)^2 */
__inline static double md_potbond(double *a, double *b,
    double r0, double k, double *fa, double *fb)
{
  double dx[D], r, dr, amp;

  r = vnorm( vdiff(dx, a, b) );
  dr = r - r0;
  if ( fa != NULL ) {
    amp = 2 * k * dr / r;
    vsinc(fa, dx, -amp);
    vsinc(fb, dx,  amp);
  }
  return k * dr * dr;
}

/* harmonic angle k (ang - ang0)^2 */
__inline static double md_potang(double *a, double *b, double *c,
    double ang0, double k, double *fa, double *fb, double *fc)
{
  double dang, amp, ga[D], gb[D], gc[D];

  dang = vang(a, b, c, ga, gb, gc) - ang0;
  if ( fa != NULL ) {
    amp = -2 * k * dang;
    vsinc(fa, ga, amp);
    vsinc(fb, gb, amp);
    vsinc(fc, gc, amp);
  }
  return k * dang * dang;
}



#if D == 3
/* 1-3 dihedral: k1 * (1 - cos(dang)) + k3 * (1 - cos(3*dang)) */
__inline static double md_potdih13(double *a, double *b, double *c, double *d,
    double ang0, double k1, double k3,
    double *fa, double *fb, double *fc, double *fd)
{
  double dang, amp, ga[3], gb[3], gc[3], gd[3], u;

  if ( fa != NULL ) {
    dang = vdih(a, b, c, d, ga, gb, gc, gd) - ang0;
    amp  = -k1 * sin(dang);
    amp += -3 * k3 * sin(3*dang);
    vsinc(fa, ga, amp);
    vsinc(fb, gb, amp);
    vsinc(fc, gc, amp);
    vsinc(fd, gd, amp);
  } else {
    dang = vdih(a, b, c, d, NULL, NULL, NULL, NULL) - ang0;
  }
  u  = k1 * (1 - cos(dang));
  u += k3 * (1 - cos(3 * dang));
  return u;
}

#endif


#endif /* MDUTIL_H__ */

