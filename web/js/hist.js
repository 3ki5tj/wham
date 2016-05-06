/* one-dimensional histograms
 * Note that the array structure is different from the C version */



"use strict";



var HIST_VERBOSE    = 0x0001;
var HIST_ADDAHALF   = 0x0010;
var HIST_NOZEROES   = 0x0020;
var HIST_KEEPLEFT   = 0x0040;
var HIST_KEEPRIGHT  = 0x0080;
var HIST_KEEPLEFT2  = 0x0040;
var HIST_KEEPRIGHT2 = 0x0080;
var HIST_KEEPEDGE   = (HIST_KEEPLEFT | HIST_KEEPRIGHT | HIST_KEEPLEFT2 | HIST_KEEPRIGHT2);
var HIST_KEEPHIST   = 0x0100;
var HIST_OVERALL    = 0x0200;
var HIST_INT        = 0x0400;
var HIST_ADDITION   = 0x1000;



function Hist(rows, xmin, xmax, dx)
{
  this.rows = rows;
  this.xmin = xmin;
  this.dx = dx;
  this.n = Math.floor((xmax - xmin)/dx + 0.99999999);
  var i, j, rows = this.rows, n = this.n;
  this.arr = newarr(this.rows);
  for ( i = 0; i < rows; i++ ) {
    this.arr[i] = newarr(n);
    for ( j = 0; j < n; j++ ) {
      this.arr[i][j] = 0;
    }
  }
}


Hist.prototype.clear = function()
{
  var rows = this.rows, n = this.n, i, j;
  for ( i = 0; i < rows; i++ ) {
    for ( j = 0; j < n; j++ ) {
      this.arr[i][j] = 0;
    }
  }
}



/* compute sum, average and variance of the histogram h */
function hist_getsums(h, n, xmin, dx)
{
  var i, x, w;
  var s = newarr(3);

  s[0] = s[1] = s[2] = 0.;
  for ( i = 0; i < n; i++ ) {
    x = xmin + (i + .5) * dx;
    w = h[i];
    s[0]  += w;
    s[1]  += w * x;
    s[2]  += w * x * x;
  }
  if ( s[0] > 0 ) {
    s[1] /= s[0];
    s[2] = s[2] - s[1] * s[1] * s[0];
    if ( s[0] > 1 ) {
      s[2] /= s[0] - 1;
    }
  }
  return s;
}



/* compute sum, average and variance of the rth histogram */
Hist.prototype.getave = function(r)
{
  return hist_getsums(this.arr[r], this.n, this.xmin, this.dx);
}



/* add x[r] of weight w into the rth row of the histogram, h, r = 0..rows-1
 * return the number of successful rows */
function histadd(xarr, w, h, rows, n, xmin, dx, flags)
{
  var r, ix, good = 0, verbose = flags & HIST_VERBOSE;
  var x;

  for (r = 0; r < rows; r++) {
    if ( (x = xarr[r]) < xmin ) {
      if (verbose)
        console.log("histadd underflows ", r, ": ", x, " < ", xmin);
      continue;
    }
    ix = Math.floor((x - xmin)/dx);
    if ( ix >= n ) {
      if (verbose)
        console.log("histadd overflows ", r, ": ", x, " > ", xmin + dx * n);
      continue;
    }
    h[r][ix] += w;
    good++;
  }
  return good;
}



/* add x[r] of weight w into the rth row of the histogram */
Hist.prototype.add1 = function(r, x, w, flags)
{
  var ix, n = this.n, verbose = flags & HIST_VERBOSE;
  var xmin = this.xmin, dx = this.dx;

  if ( r >= this.rows || r < 0 ) {
    console.log("bad row index ", r);
    return -1;
  }
  if ( x < xmin ) {
    if ( verbose )
      console.log("histadd underflows ", r, ": ", x, " < ", xmin);
    return -1;
  }
  ix = Math.floor((x - xmin)/dx);
  if ( ix >= n ) {
    if ( verbose )
      console.log("histadd overflows ", r, ": ", x, " > ", xmin+dx*n);
    return -1;
  }
  this.arr[r][ix] += w;
  return 0;
}



/* add x[r] of weight w into the rth histogram, r = 0..rows-1
 * return the number of successes */
Hist.prototype.add = function(x, w, flags)
{
  var r, good = 0;

  for ( r = 0; r < this.rows; r++ )
    good += (this.add1(r, x[r], w, flags) == 0);
  return good;
}




