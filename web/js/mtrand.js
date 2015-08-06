/* Mersenne Twister was developed by Makoto Matsumoto and Takuji Nishimura */



"use strict";



var MT_N = 624;
var MT_M = 397;
var MT_UMASK = 0x80000000;
var MT_LMASK = 0x7fffffff;

var mtonce = 0;
var mtidx = MT_N;
var mtarr = new Array(MT_N);
mtarr[0] = 5489; // the seed



/* scramble the random number state */
function mtscramble(seed)
{
  mtarr[0] = seed * 314159265 + 271828183;
  for (var k = 1; k < MT_N; k++) { // the final mask is for 64-bit machines
    mtarr[k] = 1812433253 * (mtarr[k - 1] ^ (mtarr[k - 1] >>> 30)) + k;
    // mr->arr[k] = (mr->arr[k] + seed) * 22695477ul + 1ul;
    mtarr[k] = (mtarr[k] + seed) * 314159265 + 1;
  }
  mtidx = MT_N; // request an update
  mtonce = 1; // scrambled
}



/* return an unsigned random number */
function mtrand()
{
  var mag01 = [0, 0x9908b0df]; // MATRIX_A
  var x, k;

  if ( !mtonce ) mtscramble(new Date().getTime());

  if (mtidx >= MT_N) { // generate MT_N words at one time
    for (k = 0; k < MT_N - MT_M; k++) {
      x = (mtarr[k] & MT_UMASK) | (mtarr[k+1] & MT_LMASK);
      mtarr[k] = mtarr[k+MT_M] ^ (x >>> 1) ^ mag01[x & 0x1];
    }
    for (; k < MT_N-1; k++) {
      x = (mtarr[k] & MT_UMASK) | (mtarr[k+1] & MT_LMASK);
      mtarr[k] = mtarr[k+(MT_M-MT_N)] ^ (x >>> 1) ^ mag01[x & 0x1];
    }
    x = (mtarr[MT_N-1] & MT_UMASK) | (mtarr[0] & MT_LMASK);
    mtarr[MT_N-1] = mtarr[MT_M-1] ^ (x >>> 1) ^ mag01[x & 0x1];
    mtidx = 0;
  }
  x = mtarr[ mtidx++ ];
  // tempering
  x ^= (x >>> 11);
  x ^= (x <<  7) & 0x9d2c5680;
  x ^= (x << 15) & 0xefc60000;
  x ^= (x >>> 18);
  return x >>> 0;
}



function rand01()
{
  return mtrand() / 4294967296.0;
}



/* Gaussian distribution with zero mean and unit variance
 * using the ratio method */
function randgaus()
{
  var x, y, u, v, q;
  do {
    u = 1 - rand01();
    v = 1.7156*(rand01() - .5);  // >= 2*sqrt(2/e)
    x = u - 0.449871;
    y = Math.abs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*Math.log(u));
  return v/u;
}



/* return a random number that satisfies the gamma distribution
 * p(x) = x^(k - 1) exp(-x) / Gamma(k) */
function randgam(k)
{
  var lt1 = false;
  var a, b, x, v, u;

  if ( k <= 0 ) return 0;
  if ( k < 1 ) {
    lt1 = true;
    k += 1;
  }
  a = k - 1./3;
  b = 1./3/Math.sqrt(a);

  for ( ; ; ) {
    do {
      x = randgaus();
      v = 1 + b * x;
    } while ( v <= 0 );
    v *= v * v;
    x *= x;
    u = rand01();
    if ( u <= 1 - 0.331 * x * x ) break;
    u = Math.log(u);
    if ( u <= 0.5 * x + a * (1 - v + Math.log(v)) ) break;
  }

  x = a * v;
  if ( lt1 ) x *= Math.pow(1 - rand01(), 1./(k - 1));
  return x;
}



/* return a random number that satisfies the chi-squared distribution,
 * which is the sum of the squares k Gaussian random numbers */
function randchisqr(k)
{
  return 2*randgam(k*.5);
}



/* a randomly oriented unit vector */
function randdir()
{
  var a, b, sq, s;

  do {
    a = 2 * rand01() - 1;
    b = 2 * rand01() - 1;
    sq = a * a + b * b;
  } while ( sq >= 1 );
  s = 2. * Math.sqrt(1 - sq);
  return [a * s, b * s, 1 - 2 * sq];
}

