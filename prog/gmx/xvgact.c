/* compute the autocorrelation time */
#include "../whammodel.h"
#include "xvg.h"



/* compute the autocorrelation time */
static int calcact(const char *fn, model_t *m)
{
  xvg_t *xvg;
  int k;
  double act[10];

  if ( (xvg = xvg_load(fn)) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  xvg_act(xvg, act, m->actmax, m->acmin, m->fnac);

  for ( k = 0; k < xvg->m; k++ ) {
    printf("%d: autocorrelation time %g\n", k, act[k]);
  }
  return 0;
}



int main(int argc, char **argv)
{
  model_t m[1];

  model_default(m);
  model_doargs(m, argc, argv);

  calcact(m->fninp, m);
  return 0;
}
