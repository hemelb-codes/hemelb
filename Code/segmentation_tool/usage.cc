#include "usage.h"

void Usage::visUsage (char *progname)
{
  printf ("Usage: %s input file, output config, pars, coords and checkpoint names,\n", progname);
  printf ("and checkpoint flag\n");
}
