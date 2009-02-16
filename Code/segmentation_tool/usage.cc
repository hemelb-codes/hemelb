#include "usage.h"


void visUsage (char *progname)
{
  printf ("Usage: %s input path, output config, pars and checkpoint names\n", progname);
  printf ("slice and pixel size (mm)\n");
  printf ("    or\n");
  printf ("checkpoint file name, output config and pars file names\n");
}
