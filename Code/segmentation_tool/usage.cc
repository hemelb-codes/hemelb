#include "usage.h"


void visUsage (char *progname)
{
#ifndef MESH
  printf ("Usage: %s input path, output config, pars and checkpoint names,\n", progname);
  printf ("slice and pixel size (mm), and checkpoint flag\n");
#else
  printf ("Usage: %s input file, output config, pars and checkpoint names,\n", progname);
  printf ("and checkpoint flag\n");
#endif
}
