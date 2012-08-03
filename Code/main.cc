// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
