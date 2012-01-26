#include "configuration/CommandLine.h"
#include "SimulationMaster.h"

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
