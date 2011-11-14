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



int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

  hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc,argv);
  SimulationMaster lMaster = SimulationMaster(options);

  lMaster.RunSimulation();

  return (0);
}

