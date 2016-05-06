#ifndef LJEOS_H__
#define LJEOS_H__



#include <stdio.h>
#include <math.h>



enum { LJEOS_MBWRJZG = 0, LJEOS_MBWRKN = 1, LJEOS_PVEhBHKN = 2};

#define ljeos3d_MBWRJZG(rho, T, P, Fex, muex) \
  ljeos3d_MBWR(rho, T, P, Fex, muex, LJEOS_MBWRJZG, NULL)

#define ljeos3d_MBWRKN(rho, T, P, Fex, muex) \
  ljeos3d_MBWR(rho, T, P, Fex, muex, LJEOS_MBWRKN, NULL)



/* compute reference thermal dynamics variables
   using the modified Benedic-Webb-Rubin (MBWR) equation of states
   return the average potential energy
   *P:  pressure
   *Fex: Helmholtz free energy (potential part)
   *muex: Gibbs free energy (potential part) */
__inline double ljeos3d_MBWR(double rho, double T, double *P,
    double *Fex, double *muex, int eos, const double *x)
{
  /* Reference:
   * J. Karl Johnson, John A. Zollweg, and Keith E. Gubbins (JZG)
   * The Lennard-Jones equation of states revisited,
   * Molecular Physics (1993) Vol. 78, No 3, 591-618
   * Table 10 */
  const double ljparam_JZG1993[] = {0,
     0.8623085097507421,
     2.976218765822098,
    -8.402230115796038,
     0.1054136629203555,
    -0.8564583828174598,     /* 5 */
     1.582759470107601,
     0.7639421948305453,
     1.753173414312048,
     2.798291772190376e+03,
    -4.8394220260857657e-02, /* 10 */
     0.9963265197721935,
    -3.698000291272493e+01,
     2.084012299434647e+01,
     8.305402124717285e+01,
    -9.574799715203068e+02,  /* 15 */
    -1.477746229234994e+02,
     6.398607852471505e+01,
     1.603993673294834e+01,
     6.805916615864377e+01,
    -2.791293578795945e+03,  /* 20 */
    -6.245128304568454,
    -8.116836104958410e+03,
     1.488735559561229e+01,
    -1.059346754655084e+04,
    -1.131607632802822e+02,  /* 25 */
    -8.867771540418822e+03,
    -3.986982844450543e+01,
    -4.689270299917261e+03,
     2.593535277438717e+02,
    -2.694523589434903e+03,  /* 30 */
    -7.218487631550215e+02,
     1.721802063863269e+02,
     3.0 /* gamma, Cf. JZG 1993, Table 7 caption */
  };

  /* Reference:
   * Jiri Kolafa and Ivo Nezbeda
   * The Lennard-Jones fluid: An accurate analytic
   *   and theoretically-based equation of state
   * Fluid Phase Equilibria (1994) Vol. 100, 1-34
   * TABLE 5
   * regressed from data with T <= 6
   * */
  const double ljparam_KN1994[] = {0,
        0.86230851,
        2.97621877,
       -8.40223012,
        0.10541366,
       -0.85645838, /* 5 */
        1.39945300,
       -0.20682219,
        2.66555449,
     1205.90355811,
        0.24414200, /* 10 */
        6.17927577,
      -41.33848427,
       15.14482295,
       88.90243729,
    -2425.74868591, /* 15 */
     -148.52651854,
       68.73779789,
     2698.26346845,
    -1216.87158315,
    -1199.67930914, /* 20 */
       -7.28265251,
    -4942.58001124,
       24.87520514,
    -6246.96241113,
     -235.12327760, /* 25 */
    -7241.61133138,
     -111.27706706,
    -2800.52326352,
     1109.71518240,
     1455.47321956, /* 30 */
    -2577.25311109,
      476.67051504,
        4.52000000  /* gamma */
  };

  double gam;
  double a[8], b[6], c[8], d[6], G[6];
  double F, rhop, rho2 = rho*rho, Pa = 0., Pb = 0., Pt, U, Aex;
  int i;

  /* default parameters */
  if ( eos == LJEOS_MBWRJZG ) x = ljparam_JZG1993;
  else if (eos == LJEOS_MBWRKN ) x = ljparam_KN1994;

  gam = x[33];
  /* Cf. JZG 1993, Table 5 */
  a[0] = x[1]*T + x[2]*sqrt(T) + x[3] + x[4]/T + x[5]/(T*T);
  a[1] = x[6]*T + x[7] + x[8]/T + x[9]/(T*T);
  a[2] = x[10]*T + x[11] + x[12]/T;
  a[3] = x[13];
  a[4] = x[14]/T + x[15]/(T*T);
  a[5] = x[16]/T;
  a[6] = x[17]/T + x[18]/(T*T);
  a[7] = x[19]/(T*T);
  /* Cf. JZG 1993, Table 6 */
  b[0] = (x[20] + x[21]/T)/(T*T);
  b[1] = (x[22] + x[23]/(T*T))/(T*T);
  b[2] = (x[24] + x[25]/T)/(T*T);
  b[3] = (x[26] + x[27]/(T*T))/(T*T);
  b[4] = (x[28] + x[29]/T)/(T*T);
  b[5] = (x[30] + x[31]/T + x[32]/(T*T))/(T*T);
  /* Cf. JZG 1993, Table 8 */
  c[0] = x[2]*sqrt(T)/2 + x[3] + 2*x[4]/T + 3*x[5]/(T*T);
  c[1] = x[7] + 2*x[8]/T + 3*x[9]/(T*T);
  c[2] = x[11] + 2*x[12]/T;
  c[3] = x[13];
  c[4] = 2*x[14]/T + 3*x[15]/(T*T);
  c[5] = 2*x[16]/T;
  c[6] = 2*x[17]/T + 3*x[18]/(T*T);
  c[7] = 3*x[19]/(T*T);
  /* Cf. JZG 1993, Table 9 */
  d[0] = (3*x[20] + 4*x[21]/T)/(T*T);
  d[1] = (3*x[22] + 5*x[23]/(T*T))/(T*T);
  d[2] = (3*x[24] + 4*x[25]/T)/(T*T);
  d[3] = (3*x[26] + 5*x[27]/(T*T))/(T*T);
  d[4] = (3*x[28] + 4*x[29]/T)/(T*T);
  d[5] = (3*x[30] + 4*x[31]/T + 5*x[32]/(T*T))/(T*T);

  /* Cf. JZG 1993, Table 7 */
  F = exp(-gam*rho*rho);
  G[0] = (1 - F)/(2*gam);
  for (rhop = 1, i = 1; i < 6; i++) {
    rhop *= rho*rho;
    G[i] = -(F*rhop - 2*i*G[i-1])/(2*gam);
  }

  Aex = 0;
  Pa = Pb = 0;
  for (U = 0, i = 7; i >= 0; i--) {
    /* Sum {i = 0 to 7} c[i] rho^{i+1} / (i + 1)
     * Cf. JZG 1993, Eq. (9), first term */
    U = rho * (c[i]/(i+1) + U);
    /* Sum {i = 0 to 7} a[i] rho^{i+1} / (i + 1)
     * Cf. JZG 1993, Eq. (5), first term */
    Aex = rho * (a[i]/(i+1) + Aex);
    /* rho * Sum {i = 0 to 7} a[i] rho^{i+1}
     * the leading factor rho is to be multiplied later
     * Cf. JZG 1993, Eq. (7), second term */
    Pa  = rho * (a[i] + Pa);
  }

  for (i = 5; i >= 0; i--) {
    /* Sum {i = 0 to 5} d[i] G[i]
     * Cf. JZG 1993, Eq. (9), second term */
    U += d[i]*G[i];
    /* Sum {i = 0 to 5} b[i] G[i]
     * Cf. JZG 1993, Eq. (5), second term */
    Aex += b[i]*G[i];
    /* (rho*F) * Sum {i = 0 to 5} b[i] rho^{2*(i+1)}
     * the leading factor (rho*F) is to be multiplied later
     * Cf. JZG 1993, Eq. (7), last term */
    Pb = rho2*(b[i] + Pb);
  }
  Pt = rho*(T + Pa + F*Pb); /* Cf. JZG 1993, Eq. (7) */
  if (Fex) *Fex = Aex;
  if (P) *P = Pt;
  if (muex) *muex = Aex + Pt/rho - T;
  return U;
}



__inline double ljeos3d_bAhs(double eta)
{
  double e1 = 1 - eta;
  return log(e1)*5/3 + eta*(34 - 33*eta + 4*eta*eta)/(6*e1*e1);
}



__inline double ljeos3d_zhs(double eta)
{
  double e1 = 1 - eta;
  return (1 + eta*(1 + eta*(1 - eta*(1 + eta)*2/3)))/(e1*e1*e1);
}



/* equivalent hard-sphere diameter of the hybrid Barker-Henderson theory
 * defined in Eq. (17)
 * d = Int {0 to 2^(1/6)} (1 - exp(-beta*u)) dr
 * Parameters are given by Table 2 with the functional form given by (29)
 *
 * To verify it WolframAlpha or Mathematica
 * Integrate[(1-Exp[-(4/x^12-4/x^6+1)/T]), {x, 0, 2.0^(1.0/6)}]
 * */
__inline double ljeos3d_dhBH(double T, double *dfdbeta)
{
  const double Ci[3] = {1.080142248, -0.076383859, 0.011117524};
  const double C1 = 0.000693129, Cln = -0.063920968;
  double x;

  x = sqrt(T);
  *dfdbeta = Ci[2] + 0.5*Ci[1]*x - 0.5*C1*x*T - T*Cln;
  return Ci[2]/T + Ci[1]/x + Ci[0] + C1*x + Cln*log(T);
}



/* The residual second virial coefficient B2_LJ - B2_hs
 * from the hybrid Barker-Henderson theory
 * Parameters are given by Table 2 with the functional form given by (29) */
__inline double ljeos3d_dB2hBH(double T, double *dfdbeta)
{
  /* from Table 2 */
  const double Ci[8] = {
    0.02459877,
    0,
   -7.02181962,
    2.90616279,
   -4.13749995,
    0.87361369,
    0.43102052,
   -0.58544978};
  double f = 0, invx = 1/sqrt(T);
  int i;

  *dfdbeta = 0;
  for ( i = 7; i > 0; i-- ) {
    f = (f + Ci[i]) * invx;
    *dfdbeta = (*dfdbeta - Ci[i]*i/2) * invx;
  }
  *dfdbeta *= -T;
  f += Ci[0];
  return f;
}



/* Reference:
 * Jiri Kolafa and Ivo Nezbeda
 * The Lennard-Jones fluid: An accurate analytic
 *   and theoretically-based equation of state
 * Fluid Phase Equilibria (1994) Vol. 100, 1-34
 * http://www.sklogwiki.org/SklogWiki/index.php/Lennard-Jones_equation_of_state
 */
__inline double ljeos3d_PVEhBH(double rho, double T, double *P, double *A, double *mu)
{
  /* Table 3 */
  const double Cij[5][7] = {
    {0, 0, /* i = 0 */
     2.01546797,
   -28.17881636,
    28.28313847,
   -10.42402873,
    0},
    {0, 0, /* i = -1 */
   -19.58371655,
    75.62340289,
  -120.70586598,
    93.92740328,
   -27.37737354},
    {0, 0, /* i = -2 */
    29.34470520,
  -112.35356937,
   170.64908980,
  -123.06669187,
    34.42288969},
    {0, 0, 0, 0, 0, 0, 0},
    {0, 0, /* i = -4 */
   -13.37031968,
    65.38059570,
  -115.09233113,
    88.91973082,
   -25.62099890}};
  const double gam = 1.92907278;
  double eta, Ahs, AhBH, ACij, z, zhs, zCij, U, UCij;
  double dB2, ddB2, dia, ddia, grho2, xpngrho2;
  double invx = 1/sqrt(T);
  double xprho, xpT, xp;
  int i, j;

  dia = ljeos3d_dhBH(T, &ddia);
  dB2 = ljeos3d_dB2hBH(T, &ddB2);
  grho2 = gam*rho*rho;
  xpngrho2 = exp(-grho2);
  eta = M_PI*rho*dia*dia*dia/6;
  Ahs = T*ljeos3d_bAhs(eta);
  AhBH = xpngrho2*rho*T*dB2;
  ACij = zCij = UCij = 0;
  for ( i = 0; i <= 4; i++ ) {
    xpT = pow(invx, i);
    xprho = rho;
    for ( j = 2; j <= 6; j++ ) {
      xprho *= rho;
      xp = xpT * xprho;
      ACij += Cij[i][j] * xp;
      zCij += j * Cij[i][j] * xp / T;
      UCij += (1 + i*.5) * Cij[i][j] * xp;
    }
  }
  if ( A != NULL ) *A = Ahs + AhBH + ACij;
  zhs = ljeos3d_zhs(eta);
  z = zhs + rho * (1 - 2*grho2)*xpngrho2*dB2 + zCij;
  if ( P != NULL ) *P = rho * T * z;
  U = 3*(zhs - 1)*ddia/dia + rho*xpngrho2*ddB2 + UCij;
  if ( mu != NULL ) *mu = T*(z - 1) + Ahs + AhBH + ACij;
  return U;
}



/* default */
#define ljeos3d_get(rho, T, P, Fex, muex) \
  ljeos3d_getx(rho, T, P, Fex, muex, LJEOS_PVEhBHKN)

__inline double ljeos3d_getx(double rho, double T, double *P,
    double *Fex, double *muex, int eos)
{
  if ( eos == LJEOS_MBWRJZG || eos == LJEOS_MBWRKN ) {
    return ljeos3d_MBWR(rho, T, P, Fex, muex, eos, NULL);
  } else if ( eos == LJEOS_PVEhBHKN ) {
    return ljeos3d_PVEhBH(rho, T, P, Fex, muex);
  } else {
    fprintf(stderr, "Error: unknown equation of state %d\n", eos);
    return 0;
  }
}



#endif /* LJEOS_H__ */
