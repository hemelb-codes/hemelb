// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/mpi.h"
#include "net/IOCommunicator.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation parameters and performance statistics are output on
  // standard output

  // Bring up MPI
  hemelb::net::MpiEnvironment mpi(argc, argv);
  hemelb::log::Logger::Init();
  try
  {
    hemelb::net::MpiCommunicator commWorld = hemelb::net::MpiCommunicator::World();

    hemelb::net::IOCommunicator hemelbCommunicator(commWorld);
    try
    {
      // Parse command line
      hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc, argv);

      // Start the debugger (if requested)
      hemelb::debug::Debugger::Init(options.GetDebug(), argv[0], commWorld);

      // Prepare main simulation object...
      hemelb::SimulationMaster<> master(options, hemelbCommunicator);

      // ..and run it.
      master.RunSimulation();
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
}
