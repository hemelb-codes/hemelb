#include <stdio.h>
#include "usage.h"

// Prints the correct command-line usage to the terminal. Argument is the programme name.
void Usage::printUsage(char *progname) {
  printf("Usage: %s path of the input files\n", progname);
  printf("number of pulsaticle cycles, time steps per cycle, \n");
  printf("voxel size in metres, #snapshots and #images.\n");
  printf("The following files must be present in the path specified:\n");
  printf("config.dat, pars.asc and rt_pars.asc\n");
}

