#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "util/fileutils.h"
#include "util/utilityFunctions.h"
#include "lb/lb.h"
#include "log/Logger.h"

#include "SimConfig.h"
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

  // This is currently where all default command-line arguments are.
  std::string lInputFile = "config.xml";
  std::string lOutputDir = "";
  unsigned int lSnapshotsPerCycle = 10;
  unsigned int lImagesPerCycle = 10;
  unsigned int lSteeringSessionId = 1;

  // There should be an odd number of arguments since the parameters occur in pairs.
  if ( (argc % 2) == 0)
  {
    if (lMaster.IsCurrentProcTheIOProc())
    {
      PrintUsage(argv[0]);
    }
    lMaster.Abort();
  }

  // All arguments are parsed in pairs, one is a "-<paramName>" type, and one
  // is the <parametervalue>.
  for (int ii = 1; ii < argc; ii += 2)
  {
    char* lParamName = argv[ii];
    char* lParamValue = argv[ii + 1];
    if (strcmp(lParamName, "-in") == 0)
    {
      lInputFile = std::string(lParamValue);
    }
    else if (strcmp(lParamName, "-out") == 0)
    {
      lOutputDir = std::string(lParamValue);
    }
    else if (strcmp(lParamName, "-s") == 0)
    {
      char * dummy;
      lSnapshotsPerCycle = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else if (strcmp(lParamName, "-i") == 0)
    {
      char *dummy;
      lImagesPerCycle = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else if (strcmp(lParamName, "-ss") == 0)
    {
      char *dummy;
      lSteeringSessionId = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else
    {
      if (lMaster.IsCurrentProcTheIOProc())
      {
        PrintUsage(argv[0]);
      }
      lMaster.Abort();
    }
  }

  hemelb::SimConfig *lSimulationConfig = hemelb::SimConfig::Load(lInputFile.c_str());
  // Actually create the directories.

  lMaster.SetupReporting(lOutputDir,lInputFile);

  lMaster.Initialise(lSimulationConfig, lImagesPerCycle, (int) lSteeringSessionId);

  lMaster.RunSimulation(lSnapshotsPerCycle, lImagesPerCycle);

  delete lSimulationConfig;

  return (0);
}

