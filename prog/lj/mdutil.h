#ifndef MDTIL_H__
#define MDUTIL_H__



#include "mtrand.h"
#include "mat.h"



/* remove the center of mass motion */
static void rmcom(double (*x)[D], const double *m, int n)
{
  int i;
  double xc[D] = {0}, mtot = 0, wt;

  for ( i = 0; i < n; i++ ) {
    wt = m ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);
  for ( i = 0; i < n; i++ ) {
    vdec(x[i], xc);
  }
}



#if D == 2



/* annihilate the total angular momentum */
static void shiftang(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i;
  double am, r2, xc[D] = {0, 0}, xi[D];
  double mtot = 0, wt;

  /* determine the center of mass */
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);

  am = r2 = 0.0;
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vdiff(xi, x[i], xc);
    am += wt * vcross(xi, v[i]);
    r2 += wt * vsqr(x[i]);
  }

  am = -am / r2;
  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    v[i][0] += -am * xi[1];
    v[i][1] +=  am * xi[0];
  }
}



#else



/* annihilate the total angular momentum
 * solve
 *   /  m (y^2 + z^2)   -m x y          -m x y        \
 *   |  -m x y          m (x^2 + z^2)   -m y z        |  c  =  L
 *   \  -m x z          -m y z          m (x^2 + y^2) /
 * use a velocity field
 *    v' = v - c x r
 *   */
static void shiftang(double (*x)[D], double (*v)[D],
    const double *m, int n)
{
  int i;
  double xc[D] = {0, 0, 0}, xi[D], ang[D], am[D] = {0, 0, 0};
  double dv[D], mat[D][D], inv[D][D];
  double xx = 0, yy = 0, zz = 0, xy = 0, zx = 0, yz = 0;
  double mtot = 0, wt;

  /* determine the center of mass */
  for ( i = 0; i < n; i++ ) {
    wt = ( m != NULL ) ? m[i] : 1.0;
    vsinc(xc, x[i], wt);
    mtot += wt;
  }
  vsmul(xc, 1.0 / mtot);

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
  minv(inv, mat);

  /* ang is the solution of M^(-1) * L */
  ang[0] = -vdot(inv[0], am);
  ang[1] = -vdot(inv[1], am);
  ang[2] = -vdot(inv[2], am);
  for ( i = 0; i < n; i++ ) {
    vdiff(xi, x[i], xc);
    vcross(dv, ang, xi);
    vinc(v[i], dv);
  }
}



#endif /* D == 3 */



/* compute the kinetic energy */
__inline static double md_ekin(double (*v)[D], const double *m, int n)
{
  int i;
  double ek = 0;

  if ( m == NULL ) {
    for ( i = 0; i < n; i++ ) {
      ek += vsqr( v[i] );
    }
  } else {
    for ( i = 0; i < n; i++ ) {
      ek += m[i] * vsqr( v[i] );
    }
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



static void md_langevin(double (*v)[D],
    const double *m, int n, double tp, double dt)
{
  int i, k;
  double s, v0;

  s = exp(-dt);
  v0 = sqrt( tp * (1 - s * s) );
  for ( i = 0; i < n; i++ ) {
    if ( m != NULL ) {
      v0 /= sqrt( m[i] );
    }
    for ( k = 0; k < D; k++ ) {
      v[i][k] = v[i][k] * s + v0 * randgaus();
    }
  }
}



#endif /* MDUTIL_H__ */

