/* two-dimensional Ising model */



"use strict";



/* initialize an lxl Ising model */
function Ising(l)
{
  var i, n;

  this.l = l;
  this.n = n = l * l;
  this.s = newarr(n);
  for ( i = 0; i < n; i++ ) {
    this.s[i] = -1;
  }
  this.M = -n;
  this.E = -2*n;
  this.proba = newarr(10);
  this.queue = newarr(n);
  this.used = newarr(n);
}



/* set transition probability */
Ising.prototype.setproba = function(bet)
{
  var x = Math.exp(-4 * bet);
  this.proba[0] = 1.0;
  this.proba[2] = x;
  this.proba[4] = x * x;
};



/* pick a random site, count neighbors with different spins */
Ising.prototype.pick = function()
{
  var id, ix, iy, l, n, ixp, ixm, iyp, iym;

  l = this.l;
  n = this.n;
  id = Math.floor( rand01() * n );
  ix = id % l;
  iy = id - ix;
  ixp = ( ix + 1 ) % l;
  ixm = ( ix + l - 1 ) % l;
  iyp = ( iy + l ) % n;
  iym = ( iy + n - l ) % n;
  this.h = this.s[id] * (
        this.s[iy + ixp] + this.s[iy + ixm]
      + this.s[iyp + ix] + this.s[iym + ix] );
  return id;
}



/* flip site id, with (-h) is the energy before the flip */
Ising.prototype.flip = function(id)
{
  this.s[id] = -this.s[id];
  this.M += this.s[id] * 2;
  this.E += this.h * 2;
  return this.E;
}


/* compute total energy and magnetization */
Ising.prototype.em = function()
{
  var l, n, i, j, e, m;

  e = m = 0;
  l = this.l;
  n = l * l;
  for ( i = 0; i < n; i += l ) {
    for ( j = 0; j < l; j++ ) {
      var id = i + j;
      var idr = i + (j + 1) % l;
      var idu = (i + l) % n + j;
      var s = this.s[id];
      var su = this.s[idu];
      var sr = this.s[idr];
      m += s;
      e += s * (su + sr);
    }
  }
  this.M = m;
  return this.E = -e;
}



/* add spin j to the queue if s[j] is different from s
 * return the spin */
Ising.prototype.addtoqueue = function(j, s, r)
{
  var sj = this.s[j];

  if ( sj == s && !this.used[j] && rand01() < r ) {
    this.queue[ this.cnt ] = j;
    this.cnt++;
    this.used[j] = 1;
  }
  return sj;
}



/* Wolff algorithm */
Ising.prototype.wolff = function(padd)
{
  var l = this.l, n = this.n, i, ix, iy, h = 0;

  // randomly selected a seed
  var id = Math.floor ( rand01() * n );
  var s = this.s[id];
  this.cnt = 0;
  this.queue[ this.cnt++ ] = id;
  for ( i = 0; i < n; i++ ) {
    this.used[i] = 0;
  }
  this.used[id] = 1;

  // go through spins in the queue
  for ( i = 0; i < this.cnt; i++ ) {
    id = this.queue[i];
    // flip the spin to correctly compute the local field,
    // which is the total magnetization of all spins
    // surrounding the cluster.
    //
    // consider a bond id-jd, with jd being a neighbor of id
    // 1) if jd does not make it to the cluster, then it
    //    lies on the border, and it contributes
    //    s[jd] to the local field
    // 2) if s[jd] == s, and will be included in the cluster
    //    in the future, it should contribute zero to the
    //    local field.  But we let it contribute s to the
    //    local field for now.  Since jd is added to the
    //    queue, when jd is considered in this loop, or
    //    when the bond jd-id is reconsidered, it will
    //    contribute an s[id] to the local field.  But at
    //    that time, s[id] = -s due to the flip here,
    //    so the total contribution would be s + (-s) = 0.
    this.s[id] = -s;
    // add neighbors of i with the same spins
    ix = id % l;
    iy = id - ix;
    h += this.addtoqueue(iy + (ix + 1) % l,     s, padd);
    h += this.addtoqueue(iy + (ix + l - 1) % l, s, padd);
    h += this.addtoqueue((iy + l) % n + ix,     s, padd);
    h += this.addtoqueue((iy + n - l) % n + ix, s, padd);
  }

  this.E += 2 * s * h;
  this.M -= 2 * s * this.cnt;
  return 0;
}



var LN_BIG =  50.0;

/* log(exp(a) + exp(b)) */
function lnadd(a, b)
{
  var c;
  if (a < b) { c = a; a = b; b = c; } // ensure a >= b
  return ((c = a - b) > LN_BIG) ? a : a + Math.log(1 + Math.exp(-c));
}

/* log(exp(a) - exp(b)), only works for a > b */
function lndif(a, b)
{
  var c;
  return ((c = a - b) > LN_BIG) ? a : a + Math.log(1 - Math.exp(-c));
}

/* log(exp(a)+b) */
function lnaddn(a, b)
{
  return (a > LN_BIG) ? a : a + Math.log(1 + b * Math.exp(-a));
}



/* exact solution of ising model */
function is2_exact(lx, ly, beta)
{
  var lxh, n, ex, f, th, sech, bet2, bsqr, log2, x;
  var lnz, lnz1, lnz2, lnz3, lnz4, dz, ddz;
  var z21, z31, z41, za1;
  var dr1, dr2, dr3, dr4, ddr1, ddr2, ddr3, ddr4;
  var g, g0, dg, ddg, dg0;
  var xn2b, sh2b, coth2b;
  var lnch2b, lncc2b, lncl, lnsl, cd, cdsqr, lnddcl;
  var r, sgn4 = 1;
  var eav, cv;

  lxh = .5 * lx;
  n = lx * ly;
  log2 = Math.log(2.0);
  bet2 = 2. * beta;
  bsqr = beta * beta;
  xn2b = Math.exp(-bet2);
  if (lx == 2 && ly == 2) { // 2x2 system
    var lnc, lnd;
    x = 8.*beta;
    lnc = lnadd(x, -x); // c = exp(x) + exp(-x)
    lnd = lnaddn(lnc, 6.);
    lnz = lnd + log2;
    eav = -8. * Math.exp(lndif(x, -x) - lnd); // -8*sinh(8*b)/(3+cosh(8*h))
    cv = bsqr * 384. * Math.exp(lnaddn(lnc, 2./3) - 2.0*lnd); // 64*(1+3cosh(8*b))/(3+cosh(8*b))^2
    return [lnz, eav, cv];
  } else if (Math.abs(beta) < 1e-6) { /* high T approx. normal branch unstable if beta < 1e-6 */
    lnz = n * (2.*lnadd(beta, -beta) - log2);
    x = 1. + xn2b;
    eav = -2. * n * (1. - xn2b)/x;
    cv = bsqr * 8.*n*xn2b/(x*x);
    return [lnz, eav, cv]; /* +n*tanh(beta)^4 */
  }

  lnz1 = lnz2 = lnz3 = lnz4 = 0;
  dr1 = dr2 = dr3 = dr4 = 0;
  ddr1 = ddr2 = ddr3 = ddr4 = 0;
  lnch2b = lnadd(bet2, -bet2) - log2;
  coth2b = 2./(1. - xn2b*xn2b) - 1.;
  lncc2b = lnch2b + Math.log(coth2b); // ln[ cosh(2b) * coth(2b) ]
  g0 = bet2 + Math.log(2./(1. + xn2b) - 1.);
  sgn4 = (g0 >= 0) ? 1 : -1;

  sh2b = 0.5*(1./xn2b - xn2b);
  dg0 = 2. + 2./sh2b;
  x = sh2b*sh2b;
  cd = 2. - 2./x; // cl' = cd * cosh(2b)
  cdsqr = cd*cd;
  lnddcl = lnaddn(lncc2b, 2.0/(x * sh2b)) + 2.*log2; // log(cl'')

  for (r = 0; r < ly; r++) { // for odd number
    lncl = lnaddn(lncc2b, -Math.cos((2.*r + 1.) * Math.PI / ly));
    lnsl = lncl + 0.5 * Math.log(1. - Math.exp(-2.*lncl));
    g = lnadd(lncl, lnsl);
    f = lxh*g;
    lnz1 += lnadd(f, -f);
    lnz2 += lndif(f, -f);

    dg = Math.exp(lnch2b - lnsl) * cd; // g' = cl'/sl;
    ex = Math.exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    x = lxh * dg;
    dr1 += x*th;
    dr2 += x/th;

    /* g''=cl''/sl - cl' ^2 *cl/sl^3; */
    ddg = Math.exp(lnddcl - lnsl);
    ddg -= Math.exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    sech = 2.0*dg/(ex + 1.0/ex); // g' * sech(0.5*lx*g)
    ddr1 += lxh*(ddg*th + lxh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex); // g' * sech(0.5*lx*g)
    ddr2 += lxh*(ddg/th - lxh*(sech*sech));

    if (r == 0) {
      g = g0;
    } else {
      lncl = lnaddn(lncc2b, -Math.cos(2.0 * Math.PI * r / ly));
      lnsl = lncl + 0.5 * Math.log(1- Math.exp(-2*lncl));
      g = lnadd(lncl, lnsl);
    }
    f = lxh*g;
    lnz3 += lnadd(f, -f); /* log [2 cosh(f)] */
    lnz4 += (f < 0) ? lndif(-f, f) : lndif(f, -f); /* avoid neg. g0 */

    ex = Math.exp(-f);
    th = 2./(1. + ex*ex) - 1.;
    dg = (r == 0) ? dg0 : Math.exp(lnch2b - lnsl)*cd;
    dr3 += lxh*dg*th;
    dr4 += lxh*dg/th;

    if (r == 0) {
      ddg = -4 * coth2b * coth2b * Math.exp(-lnch2b);
    } else {
      ddg = Math.exp(lnddcl - lnsl);
      ddg -= Math.exp(lnch2b*2. + lncl - 3.*lnsl)*cdsqr;
    }
    sech = 2.0*dg/(ex + 1.0/ex);
    ddr3 += lxh*(ddg*th + lxh*(sech*sech));
    sech = 2.0*dg/(ex - 1.0/ex);
    ddr4 += lxh*(ddg/th - lxh*(sech*sech));
  }

  z21 = Math.exp(lnz2 - lnz1);
  z31 = Math.exp(lnz3 - lnz1);
  z41 = sgn4 * Math.exp(lnz4 - lnz1);
  za1 = 1.0 + z21 + z31 + z41;
  lnz = lnz1 + Math.log(za1);
  lnz += .5 * n * Math.log(2.*sh2b) - log2;
  dz = (dr1 + z21*dr2 + z31*dr3 + z41*dr4)/za1;
  eav = - n*coth2b - dz;
  ddr1 += dr1*dr1;
  ddr2 += dr2*dr2;
  ddr3 += dr3*dr3;
  ddr4 += dr4*dr4;
  ddz = (ddr1 + z21*ddr2 + z31*ddr3 + z41*ddr4)/za1;
  cv = bsqr * (-2.*n/(sh2b*sh2b) + ddz - dz*dz);
  return [lnz, eav, cv];
}



