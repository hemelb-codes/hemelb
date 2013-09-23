// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/mpi.h"
#include "net/NetworkTopology.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation parameters and performance statistics are output on
  // standard output

  // Bring up MPI
  hemelb::net::MpiEnvironment mpi(argc, argv);
  try
  {
    // Start the debugger (no-op if HEMELB_USE_DEBUGGER is OFF)
    hemelb::debug::Debugger::Init(argv[0]);

    hemelb::net::MpiCommunicator hemelbCommunicator(MPI_COMM_WORLD);
    hemelb::net::NetworkTopology::Instance()->Init(hemelbCommunicator);
    try
    {
      // Parse command line
      hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc, argv);

      // Prepare main simulation object...
      SimulationMaster master = SimulationMaster(options);

      // ..and run it.
      master.RunSimulation();
    }

    // Interpose this catch to print usage before propagating the error.
    catch (hemelb::configuration::CommandLine::OptionError& e)
    {
      if (hemelbCommunicator.GetRank() == 0)
      {
        hemelb::configuration::CommandLine::PrintUsage();
      }
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
