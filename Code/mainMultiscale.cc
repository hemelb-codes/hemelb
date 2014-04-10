// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

/* Specifying the use of 'real' MPWide */
#include <MPWide.h>

#include "configuration/CommandLine.h"
#include "multiscale/MultiscaleSimulationMaster.h"
#include "multiscale/mpwide/MPWideIntercommunicator.h"

int main(int argc, char *argv[])
{
  // Bring up MPI
  hemelb::net::MpiEnvironment mpi(argc, argv);
  hemelb::log::Logger::Init();
  try
  {
    hemelb::net::MpiCommunicator commWorld = hemelb::net::MpiCommunicator::World();
    // Start the debugger (no-op if HEMELB_USE_DEBUGGER is OFF)
    hemelb::debug::Debugger::Init(argv[0], commWorld);

    hemelb::net::IOCommunicator& hemelbCommunicator = hemelb::net::IOCommunicator::Init(commWorld);

    try
    {
      // Parse command line
      hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc, argv);

      // Prepare some multiscale/MPWide stuff

      // Work out the location of the input file.
      std::string inputFile = options.GetInputFile();
      int sl = inputFile.find_last_of("/");
      std::string mpwideConfigDir = inputFile.substr(0, sl + 1);

      std::cout << inputFile << " " << mpwideConfigDir << " " << sl << std::endl;

      // Create some necessary buffers.
      std::map<std::string, bool> lbOrchestration;
      lbOrchestration["boundary1_pressure"] = true;
      lbOrchestration["boundary2_pressure"] = true;
      lbOrchestration["boundary1_velocity"] = true;
      lbOrchestration["boundary2_velocity"] = true;

      std::map<std::string, double> sharedValueBuffer;

      // TODO The MPWide config file should be read from the HemeLB XML config file!
      // Create the intercommunicator
      hemelb::multiscale::MPWideIntercommunicator intercomms(hemelbCommunicator.OnIORank(),
                                                             sharedValueBuffer,
                                                             lbOrchestration,
                                                             mpwideConfigDir.append("MPWSettings.cfg"));

      //TODO: Add an IntercommunicatorImplementation?
      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Constructing MultiscaleSimulationMaster()");
      hemelb::multiscale::MultiscaleSimulationMaster<hemelb::multiscale::MPWideIntercommunicator> lMaster(options,
                                                                                                          intercomms);

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Runing simulation()");
      lMaster.RunSimulation();
    }
    // Interpose this catch to print usage before propagating the error.
    catch (hemelb::configuration::CommandLine::OptionError& e)
    {
      hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::Singleton>(hemelb::configuration::CommandLine::GetUsage());
      throw;
    }
  }
  catch (std::exception& e)
  {
    hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::OnePerCore>(e.what());
    mpi.Abort(-1);
  }
  // MPI gets finalised by MpiEnv's d'tor.
  return (0);
}
