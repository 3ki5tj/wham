/* reference values for two-dimensional Ising model */
#define IS2_LB 6
#include "is2.h"
#include <time.h>
#define IS2_MODEL
#include "../whammodel.h"



static void model_default_is2(model_t *m)
{
  model_default(m);
  m->nT = 1600;
  m->Tmin = 1.5;
  m->Tdel = 0.001;
}



int main(int argc, char **argv)
{
  model_t m[1];
  double T, lnz, lnz0 = 0;
  int iT;

  model_default_is2(m);
  model_doargs(m, argc, argv);

  for ( iT = 0; iT < m->nT; iT++ ) {
    T = m->Tmin + m->Tdel * iT;
    lnz = is2_exact(IS2_L, IS2_L, 1/T, NULL, NULL);
    if ( iT == 0 ) {
      lnz0 = lnz;
    }
    printf("%5.3f %12.7f\n", T, lnz - lnz0);
  }

  return 0;
}

