#ifndef MTRAND_H__
#define MTRAND_H__



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



/* Mersenne Twister was developed by Makoto Matsumoto and Takuji Nishimura */
#define MT_N 624
#define MT_M 397
#define MT_UMASK 0x80000000u /* most significant w-r bits */
#define MT_LMASK 0x7fffffffu /* least significant r bits */



int mtonce = 0;
int mtidx = MT_N;
unsigned long mtarr[MT_N] = {5491};



/* scramble the random number state */
__inline static void mtscramble(unsigned seed)
{
  int k;

  mtarr[0] = (seed * 314159265u + 271828183u) & 0xfffffffful;
  for (k = 1; k < MT_N; k++) { /* the final mask is for 64-bit machines */
    mtarr[k] = 1812433253ul * (mtarr[k - 1] ^ (mtarr[k - 1] >> 30)) + k;
    /* mr->arr[k] = (mr->arr[k] + seed) * 22695477ul + 1ul; */
    mtarr[k] = ((mtarr[k] + seed) * 314159265ul + 1ul) & 0xfffffffful;
  }
  mtidx = MT_N; /* request for an update */
  mtonce = 1; /* scrambled */
}



/* return an unsigned random number */
__inline static unsigned mtrand(void)
{
  static const unsigned mag01[2] = {0, 0x9908b0dfu}; /* MATRIX_A */
  unsigned x;
  int k;

  if ( !mtonce ) mtscramble(12345);

  if (mtidx >= MT_N) { /* generate MT_N words at one time */
    for (k = 0; k < MT_N - MT_M; k++) {
      x = (mtarr[k] & MT_UMASK) | (mtarr[k+1] & MT_LMASK);
      mtarr[k] = mtarr[k+MT_M] ^ (x>>1) ^ mag01[x&1u];
    }
    for (; k < MT_N-1; k++) {
      x = (mtarr[k] & MT_UMASK) | (mtarr[k+1] & MT_LMASK);
      mtarr[k] = mtarr[k+(MT_M-MT_N)] ^ (x>>1) ^ mag01[x&1u];
    }
    x = (mtarr[MT_N-1] & MT_UMASK) | (mtarr[0] & MT_LMASK);
    mtarr[MT_N-1] = mtarr[MT_M-1] ^ (x>>1) ^ mag01[x&1u];
    mtidx = 0;
  }
  x = mtarr[ mtidx++ ];
  /* tempering */
  x ^= (x >> 11);
  x ^= (x <<  7) & 0x9d2c5680u;
  x ^= (x << 15) & 0xefc60000u;
  x ^= (x >> 18);
  return x;
}



__inline static double rand01(void)
{
  return mtrand() / 4294967296.0;
}



/* Gaussian distribution with zero mean and unit variance
 * using the ratio method */
__inline static double randgaus(void)
{
  double x, y, u, v, q;
  do {
    u = 1 - rand01();
    v = 1.7156*(rand01() - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));
  return v/u;
}



/* return a random number that satisfies the gamma distribution
 * p(x) = x^(k - 1) exp(-x) / Gamma(k) */
__inline static double randgam(double k)
{
  int lt1 = 0;
  double a, b, x, v, u;

  if ( k <= 0 ) return 0;
  if ( k < 1 ) {
    lt1 = 1;
    k += 1;
  }
  a = k - 1./3;
  b = 1./3/sqrt(a);

  for ( ; ; ) {
    do {
      x = randgaus();
      v = 1 + b * x;
    } while ( v <= 0 );
    v *= v * v;
    x *= x;
    u = rand01();
    if ( u <= 1 - 0.331 * x * x ) break;
    u = log(u);
    if ( u <= 0.5 * x + a * (1 - v + log(v)) ) break;
  }

  x = a * v;
  if ( lt1 ) x *= pow(1 - rand01(), 1./(k - 1));
  return x;
}



/* return a random number that satisfies the chi-squared distribution,
 * which is the sum of the squares k Gaussian random numbers */
__inline static double randchisqr(double k)
{
  return 2*randgam(k*.5);
}



/* a randomly oriented unit vector */
__inline static double *randdir(double *v)
{
  double a, b, sq, s;

  do {
    a = 2 * rand01() - 1;
    b = 2 * rand01() - 1;
    sq = a * a + b * b;
  } while ( sq >= 1 );
  s = 2. * sqrt(1 - sq);
  v[0] = a * s;
  v[1] = b * s;
  v[2] = 1 - 2 * sq;
  return v;
}



/* save the current state to file */
__inline static int mtsave(const char *fn)
{
  FILE *fp;
  int k;

  if ( !mtonce ) return 1; /* never used */
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "MTSEED\n%d\n", mtidx);
  for (k = 0; k < MT_N; k++)
    fprintf(fp, "%lu\n", mtarr[k]);
  fclose(fp);
  return 0;
}



/* load state from `fn' */
__inline static int mtload(const char *fn, unsigned long seed)
{
  char s[64];
  int k, z, err = 1;
  FILE *fp;

  if ( (fp = fopen(fn, "r")) != NULL ) { /* try to load from file */
    if ( fgets(s, sizeof s, fp) == NULL ) {
      fprintf(stderr, "%s is empty\n", fn);
    } else if ( strncmp(s, "MTSEED", 6) != 0 ) { /* to check the first line */
      fprintf(stderr, "%s corrupted\n", fn);
    } else if ( fscanf(fp, "%d", &mtidx) != 1 ) {
      fprintf(stderr, "no index in %s\n", fn);
    } else {
      if ( !mtonce ) mtidx = MT_N; /* request updating */
      for ( z = 1, k = 0; k < MT_N; k++ ) {
        if ( fscanf(fp, "%lu", &mtarr[k]) != 1 ) break;
        if ( mtarr[k] != 0 ) z = 0; /* a non-zero number */
      }
      if ( k != MT_N ) {
        fprintf(stderr, "%s incomplete %d/%d\n", fn, k, MT_N);
      } else {
        err = z; /* clear error, if array is nonzero */
        mtonce = 1;
      }
    }
    fclose(fp);
  }

  if (err) mtscramble(seed);
  return !mtonce;
}



#endif
