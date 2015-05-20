/* compute the autocorrelation time */
#include "xvgmodel.h"
#include "xvg.h"



/* compute the autocorrelation time */
static int calcact(const char *fn, model_t *m)
{
  xvg_t *xvg;
  double act;

  if ( (xvg = xvg_load(fn)) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  act = xvg_act(xvg, m->actmax, m->acmin, m->fnac);

  printf("autocorrelation time %g\n", act);
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
