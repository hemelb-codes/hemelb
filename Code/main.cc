#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "util/fileutils.h"
#include "util/utilityFunctions.h"
#include "lb/lb.h"
#include "log/Logger.h"

#include "SimulationMaster.h"

#include "debug/Debugger.h"

void PrintUsage(char *progname)
{
  printf("-!-!-!-!-!-!-!-!-!-!-!-!");
  printf("Correct usage: %s [-<Parameter Name> <Parameter Value>]* \n", progname);
  printf("Parameter name and significance:\n");
  printf("-in \t Path to the configuration xml file (default is config.xml)\n");
  printf("-out \t Path to the output folder (default is based on input file, e.g. config_xml_results)\n");
  printf("-s \t Number of snapshots to take per cycle (default 10)\n");
  printf("-i \t Number of images to create per cycle (default is 10)\n");
  printf("-ss \t Steering session identifier (default is 1)\n");
  printf("-!-!-!-!-!-!-!-!-!-!-!-!");
}

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

  SimulationMaster lMaster = SimulationMaster(argc, argv);

  lMaster.RunSimulation();


  return (0);
}

