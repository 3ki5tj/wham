/* Tempering */



"use strict";



/* create a SimTemp object */
function SimTemp(beta)
{
  this.n = beta.length;
  this.beta = beta;
  this.hist = newarr(this.n);
  this.lnz = newarr(this.n);
  this.usum = newarr(this.n);
  this.u2sum = newarr(this.n);
  this.lnzref = newarr(this.n);
  this.uref = newarr(this.n);
  this.cvref = newarr(this.n);
}



SimTemp.prototype.jump = function(itp, ep, type)
{
  var jtp, n = this.n;

  if ( type == 1 ) {
    // choose any temperature other than `itp`
    jtp = ( itp + 1 + Math.floor(rand01() * (n - 1)) ) % n;
  } else {
    // choose one of the neighbors
    if ( rand01() > 0.5 ) {
      jtp = itp + 1;
      if ( jtp >= n ) return itp;
    } else {
      jtp = itp - 1;
      if ( jtp < 0 ) return itp;
    }
  }

  // compute the acceptance probability
  var x = (this.beta[itp] - this.beta[jtp]) * ep
        + this.lnz[itp] - this.lnz[jtp];
  if ( x > 0 ) {
    return jtp;
  } else {
    var r = rand01();
    return r < Math.exp(x) ? jtp : itp;
  }
}


SimTemp.prototype.add = function(itp, ep)
{
  this.hist[itp] += 1;
  this.usum[itp] += ep;
  this.u2sum[itp] += ep * ep;
}
