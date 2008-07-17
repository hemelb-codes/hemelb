#include <stdio.h>
#include "usage.h"

void usage(char *progname) {
  printf("Usage: %s path of the input files and minutes for benchmarking\n", progname);
  printf("if one wants to do a benchmark or\n");
  printf("number of pulsaticle cycles, time steps per cycle and\n");
  printf("voxel size in metres otherwise.\n");
  printf("The following files must be present in the path specified:\n");
  printf("config.dat, pars.asc rt_pars.asc\n");

}

