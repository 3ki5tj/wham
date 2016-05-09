/* Weighted histogram analysis method
 * this module requires `hist.h` */



"use strict";



var LOG0 = -1e9;



/*
const char *wham_methods[] = {
  "Direct",
  "MDIIS",
  "ST",
  "UI",
  "WHAM_NMETHODS"
};
*/



function WHAM(beta, hist, flags)
{
  this.beta = beta;
  this.hist = hist;

  this.nbeta = beta.length;
  var nbeta = this.nbeta;
  this.tot = newarr(nbeta);
  this.lntot = newarr(nbeta);
  this.res = newarr(nbeta);

  var n = hist.n;
  this.htot = newarr(n);
  this.lndos = newarr(n);
  this.sum = newarr(nbeta);
  this.ave = newarr(nbeta);
  this.var = newarr(nbeta);
  this.hnorm = newarr(nbeta);
  this.ave1 = newarr(nbeta);
  this.var1 = newarr(nbeta);
  this.imin = 0;
  this.imax = 0;
  this.flags = flags;

  var i, j;
  // compute the total
  for ( i = 0; i < n; i++ ) {
    this.htot[i] = 0;
  }
  for ( j = 0; j < nbeta; j++ ) {
    var x, s = 0;
    for ( i = 0; i < n; i++ ) {
      x = hist.arr[j][i];
      s += x;
      this.htot[i] += x;
    }
    this.tot[j] = s;
    this.lntot[j] = (s > 0) ? Math.log(s) : LOG0;
  }

  // determine the boundaries
  // find imin
  for ( i = 0; i < n; i++ ) {
    if ( this.htot[i] > 0 ) break;
  }
  this.imin = i;

  // find imax
  for ( i = n - 1; i >= 0; i-- ) {
    if ( this.htot[i] > 0 ) break;
  }
  this.imax = i + 1;

  console.log("tot: ", this.tot, "i: [",  this.imin, ", ", this.imax, ")");
}



/* flags */
var  WHAM_RMCOM = 0x0010; // normalize by removing the center of mass motion
var WHAM_NOEST = 0x0020; // do not estimate lnz at the beginning

/* log(exp(a) + exp(b)) */
function wham_lnadd(a, b)
{
  var c;
  if (a < b) { c = a; a = b; b = c; } // ensure a >= b
  return ((c = a - b) > 50.0) ? a : a + Math.log(1 + Math.exp(-c));
}



/* compute the energy and heat capacity from the WHAM */
WHAM.prototype.getav = function()
{
  var hist = this.hist;
  var T, b, e, ee, lne, lnz, slne, slnee, lnw;
  var de = hist.dx, emin = hist.xmin;
  var i, j, nbeta = hist.rows;

  for ( j = 0; j < nbeta; j++ ) {
    b = this.beta[j];
    T = 1 / b;
    lnz = slne = slnee = LOG0;
    for ( i = this.imin; i < this.imax; i++ ) {
      if ( this.lndos[i] <= LOG0 ) continue;
      // note: we do not add emin here for it may lead to
      // a negative energy whose logarithm is undefined
      e = (i + .5) * de;
      lne = Math.log(e);
      lnw = this.lndos[i] - b * e;
      lnz = wham_lnadd(lnz, lnw);
      slne = wham_lnadd(slne, lne + lnw);
      slnee = wham_lnadd(slnee, 2*lne + lnw);
    }
    e = Math.exp(slne - lnz);
    ee = Math.exp(slnee - lnz) - e * e;
    e += emin;
    console.log(j, b, e, ee);
  }
}



/* shift `lnz` such that `lnz[0] == 0` */
function wham_shift(lnz, nbeta)
{
  var i;

  for ( i = 1; i < nbeta; i++ ) {
    lnz[i] -= lnz[0];
  }
  lnz[0] = 0;
}



/* normalize by remove the center of mass motion */
function wham_rmcom(arr, tot, n)
{
  var i;
  var s = 0, sy = 0;

  for ( i = 0; i < n; i++ ) {
    s += tot[i];
    sy += arr[i] * tot[i];
  }
  sy /= s;
  for ( i = 0; i < n; i++ ) {
    arr[i] -= sy;
  }
}



/* normalize lnz */
function wham_normalize(lnz, nbeta, w)
{
  if ( w.flags & WHAM_RMCOM ) {
    wham_rmcom(lnz, w.tot, nbeta);
  } else {
    wham_shift(lnz, nbeta);
  }
}



/* estimate the partition function using the single histogram method */
WHAM.prototype.estimatelnz = function(lnz)
{
  var hist = this.hist;
  var i, j, n = hist.n, nbeta = hist.rows;
  var db, e, h, s, dlnz;
  var imin = this.imin, imax = this.imax;

  lnz[0] = 0;
  for ( j = 1; j < nbeta; j++ ) {
    /* estimate the free energy different between
     * the pair j - 1 and j */
    db = this.beta[j] - this.beta[j - 1];
    s = 0;
    dlnz = LOG0;
    /* Z(j-1) / Z(j) = < exp( [beta(j) - beta(j-1)] E ) >_j
     *    Sum_E h_j(E) exp( [beta(j) - beta(j-1)] E )
     * = ---------------------------------------------
     *               Sum_E h_j(E)
     **/
    for ( i = imin; i < imax; i++ ) {
      h = hist.arr[j][i];
      if ( h <= 0 ) continue;
      e = hist.xmin + (i + .5) * hist.dx;
      s += h;
      dlnz = wham_lnadd(dlnz, Math.log(h) + db * e);
    }
    lnz[j] = lnz[j - 1] + (s > 0 ? Math.log(s) - dlnz : 0);
  }

  wham_normalize(lnz, nbeta, this);
}



/* compute the partition function from the density of states */
WHAM.prototype.getlnz = function(lnz)
{
  var hist = this.hist;
  var i, j, nbeta = hist.rows;
  var e;

  var imin = this.imin, imax = this.imax;
  for ( j = 0; j < nbeta; j++ ) {
    lnz[j] = LOG0;
    for ( i = imin; i < imax; i++ ) {
      if ( this.lndos[i] <= LOG0 ) continue;
      e = hist.xmin + (i + .5) * hist.dx;
      lnz[j] = wham_lnadd(lnz[j], this.lndos[i] - this.beta[j] * e);
    }
  }
  wham_normalize(lnz, nbeta, this);
}



WHAM.prototype.step = function(lnz, res, damp)
{
  var hist = this.hist;
  var i, j, n = hist.n, nbeta = hist.rows;
  var x, lnden, e, emin = hist.xmin, de = hist.dx, err;
  var imin = this.imin, imax = this.imax;

  for ( i = imin; i < imax; i++ ) {
    if ( this.htot[i] <= 0 ) {
      this.lndos[i] = LOG0;
      continue;
    }

    e = emin + (i + .5) * de;
    lnden = LOG0;
    //        num           Sum_j h_j(i)
    // dos = ----- = ------------------------------------------
    //        den     Sum_j tot_j exp(-beta_j * e) / Z_j
    for ( j = 0; j < nbeta; j++ ) {
      lnden = wham_lnadd(lnden, this.lntot[j] - this.beta[j] * e - lnz[j]);
    }
    this.lndos[i] = Math.log(this.htot[i]) - lnden;
  }

  // shift the baseline of the density of states
  for ( x = this.lndos[imin], i = 0; i < n; i++ )
    if ( this.lndos[i] > LOG0 )
      this.lndos[i] -= x;

  // refresh the partition function, save it in `res`
  this.getlnz(res);

  for ( err = 0, i = 0; i < nbeta; i++ ) {
    res[i] -= lnz[i];
    if ( Math.abs(res[i]) > err ) {
      err = Math.abs(res[i]);
    }
  }

  if ( damp > 0 ) {
    for ( i = 0; i < nbeta; i++ ) {
      lnz[i] += damp * res[i];
    }
    wham_normalize(lnz, nbeta, this);
  }

  return err;
}



/* iteratively compute the logarithm of the density of states
 * using the weighted histogram method */
WHAM.prototype.getlndos = function(lnz,
    damp, itmin, itmax, tol, verbose)
{
  var it;
  var err, errp;

  err = errp = 1e30;
  for ( it = 0; it < itmax; it++ ) {
    err = this.step(lnz, this.res, damp);
    if ( verbose ) {
      console.log("it " + it + ", err " + errp + " -> " + err + ", tol " + tol);
    }
    if ( err < tol && it > itmin ) {
      break;
    }
    errp = err;
  }

  return "WHAM converged in " + it + " steps,<br>error " + err + "\n";
}



/* weighted histogram analysis method */
function wham(hist, beta, lnz, flags, damp,
    itmin, itmax, tol, verbose)
{
  var w = new WHAM(beta, hist, flags);
  var err;

  if ( !(flags & WHAM_NOEST) ) {
    w.estimatelnz(lnz);
  }

  err = w.getlndos(lnz, damp, itmin, itmax, tol, verbose);
  //w.getav();
  return err;
}



/* non-iteratively compute the logarithm of the density of states
 * using the statistical-temperature WHAM
 * if `gauss`, use Gaussian correlation for missing data */
WHAM.prototype.stwham_getlndos = function()
{
  var hist = this.hist;
  var i, j, n = hist.n, nbeta = hist.rows;
  var imin = this.imin, imax = this.imax;
  var x, y, stbeta, de = hist.dx;

  for ( i = imin; i < imax; i++ ) {
    if ( this.htot[i] <= 0 ) continue;

    stbeta = 0;
    for ( j = 0; j < nbeta; j++ ) {
      /* compute the statistical temperature from copy j
       *            sum_j [(n_j)'(E) + n_j(E) beta_j]
       * beta(E) = -----------------------------------
       *                     sum_k n_k(E)
       * y is the numerator */
      stbeta += hist.arr[j][i] * this.beta[j];
    }

    // lndos currently holds the second part of
    // the statistical temperature
    this.lndos[i] = stbeta / this.htot[i];
  }

  for ( i = imin; i < imax; i++ ) {
    var il, ir;

    if ( this.htot[i] > 0 ) continue;

    /* find the closest nonempty bin */
    for ( il = i - 1; il >= 0; il-- ) {
      if ( this.htot[il] > 0 ) break;
    }
    for ( ir = i + 1; ir < n; ir++ ) {
      if ( this.htot[ir] > 0 ) break;
    }
    if ( il >= 0 && i - il <= ir - i ) {
      stbeta = this.lndos[il];
    } else if ( ir < n ) {
      stbeta = this.lndos[ir];
    } else {
      stbeta = 0;
    }
    this.lndos[i] = stbeta;
  }

  // integrate the second part of the statistical temperature
  // to get the density of states
  x = y = 0;
  for ( i = imin; i < imax; i++ ) {
    stbeta = this.lndos[i];
    this.lndos[i] = y + (x + stbeta) / 2 * de;
    x = stbeta; // previous temperature
    y = this.lndos[i]; // previous lndos
  }

  // add the first part
  for ( i = imin; i < imax; i++ ) {
    if ( this.htot[i] > 0 ) {
      this.lndos[i] += Math.log( this.htot[i] );
    } else {
      this.lndos[i] = LOG0;
    }
  }

  for ( i = 0; i < imin; i++ ) {
    this.lndos[i] = LOG0;
  }
  for ( i = imax; i < n; i++ ) {
    this.lndos[i] = LOG0;
  }

  return "ST-WHAM completed";
}



/* statistical temperature weighted histogram analysis method */
function stwham(hist, beta, lnz, flags)
{
  var w = new WHAM(beta, hist, flags);

  var err = w.stwham_getlndos();
  w.getlnz(lnz);
  //w.getav();
  console.log("ST-WHAM", w.tot, lnz);
  return err;
}



/* compute histogram averages */
WHAM.prototype.gethave = function()
{
  var hist = this.hist;
  var i, j, n = hist.n, nbeta = hist.rows;
  var tot, sx, sxx, x, y, de = hist.dx;
  var imin = this.imin, imax = this.imax;

  // compute averages and variance
  for ( j = 0; j < nbeta; j++ ) {
    tot = sx = sxx = 0;
    for ( i = imin; i < imax; i++ ) {
      y = hist.arr[j][i];
      x = hist.xmin + ( i + 0.5 ) * de;
      tot += y;
      sx += x * y;
      sxx += x * x * y;
    }
    sx /= tot;
    sxx = sxx / tot - sx * sx;
    if ( tot > 1 ) {
      sxx *= tot / (tot - 1);
    }
    this.sum[j] = tot;
    this.ave[j] = sx;
    this.var[j] = sxx;
  }
}



/* non-iteratively compute the logarithm of the density of states
 * using umbrella integration */
WHAM.prototype.umbint_getlndos = function()
{
  var hist = this.hist;
  var i, j, imin, imax, n = hist.n, nbeta = hist.rows;
  var x, y, tot, stbeta, de = hist.dx;

  /* compute averages and variance */
  this.gethave();
  for ( j = 0; j < nbeta; j++ ) {
    this.hnorm[j] = this.sum[j] * de / Math.sqrt(2 * Math.PI * this.var[j]);
  }

  imin = 0; // this.imin;
  imax = n; // this.imax;
  for ( i = imin; i < imax; i++ ) {
    var ei = hist.xmin + (i + 0.5) * de;

    tot = 0;
    stbeta = 0;
    for ( j = 0; j < nbeta; j++ ) {
      // use Gaussian approximation */
      var De = this.ave[j] - ei;
      x = this.hnorm[j] * Math.exp( -0.5 * De * De / this.var[j] );
      tot += x;
      stbeta += x * (this.beta[j] + De / this.var[j]);
    }

    // lndos currently holds the statistical temperature
    this.lndos[i] = (tot > 0) ? stbeta / tot : 0;
  }

  /* integrate the statistical temperature
   * to get the density of states */
  x = y = 0;
  for ( i = imin; i < imax; i++ ) {
    stbeta = this.lndos[i];
    this.lndos[i] = y + 0.5 * (x + stbeta) * de;
    x = stbeta; // previous temperature
    y = this.lndos[i]; // previous lndos
  }

  for ( i = 0; i < imin; i++ ) {
    this.lndos[i] = LOG0;
  }
  for ( i = imax; i < n; i++ ) {
    this.lndos[i] = LOG0;
  }

  return "Umbrella integration completed";
}



/* umbrella integration */
function umbint(hist, beta, lnz, flags)
{
  var w = new WHAM(beta, hist, flags);

  var err = w.umbint_getlndos();
  w.getlnz(lnz);
  //w.getav();
  console.log("UI", w.tot, lnz);
  return err;
}



function wham_getres(w, lnz, res)
{
  return w.step(lnz, res, 0);
}



function wham_mdiis(hist, beta, lnz, flags, lnzref,
    nbases, damp, update_method, threshold,
    itmin, itmax, tol, verbose)
{
  var w = new WHAM(beta, hist, flags);
  var err;

  if ( !(flags & WHAM_NOEST) ) {
    w.estimatelnz(lnz);
  }
  if ( lnzref ) {
    wham_normalize(lnzref, hist.rows, w);
  }

  err = iter_mdiis(lnz, hist.rows,
      wham_getres, wham_normalize, w,
      nbases, damp, update_method, threshold,
      itmin, itmax, tol, lnzref, verbose);

  //w.getav();
  return err;
}



/* convenience wrapper of WHAM */
function whamx(hist, beta, lnz, flags, lnzref,
    damp, nbases, update_method, threshold,
    itmin, itmax, tol, verbose, method)
{
  if ( method == "WHAM_DIRECT" ) {
    return wham(hist, beta, lnz, flags,
        damp, itmin, itmax, tol, verbose);
  } else if ( method == "WHAM_ST" ) {
    return stwham(hist, beta, lnz, flags);
  } else if ( method == "WHAM_UI" ) {
    return umbint(hist, beta, lnz, flags);
  } else if ( method == "WHAM_MDIIS" ) {
    return wham_mdiis(hist, beta, lnz, flags, lnzref,
        nbases, damp, update_method, threshold,
        itmin, itmax, tol, verbose);
  }

  return 0;
}



