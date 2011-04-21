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

  unsigned long lLastForwardSlash = lInputFile.rfind('/');
  if (lOutputDir.length() == 0)
  {
    lOutputDir = ( (lLastForwardSlash == std::string::npos)
      ? "./"
      : lInputFile.substr(0, lLastForwardSlash)) + "results";
  }

  double total_time = hemelb::util::myClock();

  FILE *timings_ptr = NULL;
  std::string image_directory = lOutputDir + "/Images/";
  std::string snapshot_directory = lOutputDir + "/Snapshots/";

  // Actually create the directories.


  if (lMaster.IsCurrentProcTheIOProc())
  {
    if (hemelb::util::DoesDirectoryExist(lOutputDir.c_str()))
    {
      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("\nOutput directory \"%s\" already exists. Exiting.",
                                                                          lOutputDir.c_str());
      lMaster.Abort();
    }

    hemelb::util::MakeDirAllRXW(lOutputDir);
    hemelb::util::MakeDirAllRXW(image_directory);
    hemelb::util::MakeDirAllRXW(snapshot_directory);

    // Save the computed config out to disk in the output directory so we have
    // a record of the total state used.
    std::string lFileNameComponent = std::string( (lLastForwardSlash == std::string::npos)
      ? lInputFile
      : lInputFile.substr(lLastForwardSlash));
    lSimulationConfig->Save(lOutputDir + "/" + lFileNameComponent);

    char timings_name[256];
    char procs_string[256];

    sprintf(procs_string, "%i", lMaster.GetProcessorCount());
    strcpy(timings_name, lOutputDir.c_str());
    strcat(timings_name, "/timings");
    strcat(timings_name, procs_string);
    strcat(timings_name, ".asc");
    timings_ptr = fopen(timings_name, "w");
    fprintf(timings_ptr, "***********************************************************\n");
    fprintf(timings_ptr, "Opening config file:\n %s\n", lInputFile.c_str());
  }

  lMaster.Initialise(lSimulationConfig, lImagesPerCycle, (int) lSteeringSessionId, timings_ptr);

  lMaster.RunSimulation(lSimulationConfig,
                        total_time,
                        image_directory,
                        snapshot_directory,
                        lSnapshotsPerCycle,
                        lImagesPerCycle);

  if (lMaster.IsCurrentProcTheIOProc())
  {
    fclose(timings_ptr);
  }

  delete lSimulationConfig;

  return (0);
}

