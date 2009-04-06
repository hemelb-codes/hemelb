#include <stdio.h>
#include "usage.h"

void usage(char *progname) {
  printf("Usage: %s path of the input files\n", progname);
  printf("number of time steps per cycle, voxel size in metres and\n");
  printf("minutes for benchmarking if one wants to do a benchmark or\n");
  printf("number of pulsaticle cycles, time steps per cycle, \n");
  printf("voxel size in metres, #snapshots and #images otherwise.\n");
  printf("The following files must be present in the path specified:\n");
  printf("config.dat, pars.asc and rt_pars.asc\n");

}

