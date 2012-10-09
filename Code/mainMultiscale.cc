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
#include "multiscale/MultiscaleSimulationMaster.h"
#include "resources/Resource.h"
#include "multiscale/Intercommunicator.h"

#include "MPWide.h" /* Specifying the use of 'real' MPWide */
#include "multiscale/mpwide/MPWideIntercommunicator.h"

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

  std::map<std::string,double> *sharedValueBuffer;
  std::map<std::string,bool> *lbOrchestration;

  hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc,argv);

  string test = options.GetInputFile();
  string s("/");
  int sl = test.find_last_of(s);
  string mpw_config_dir = test.substr(0,sl+1);

  std::cout << test << " " << mpw_config_dir << " " << sl << endl;
  
  hemelb::multiscale::MultiscaleSimulationMaster<hemelb::multiscale::MPWideIntercommunicator> *lMaster;

  sharedValueBuffer = new std::map<std::string, double>();

  lbOrchestration=new std::map<std::string,bool>();
  std::map<std::string,bool> &rorchestrationLB=*lbOrchestration;
  rorchestrationLB["boundary1_pressure"] = true;
  rorchestrationLB["boundary2_pressure"] = true;
  rorchestrationLB["boundary1_velocity"] = true;
  rorchestrationLB["boundary2_velocity"] = true;

  hemelb::multiscale::mpwide::mpwide_config_file = mpw_config_dir.append("MPWSettings.cfg");
  hemelb::multiscale::MPWideIntercommunicator intercomms(*sharedValueBuffer,*lbOrchestration);
  //TODO: Add an IntercommunicatorImplementation?

  lMaster = new hemelb::multiscale::MultiscaleSimulationMaster<hemelb::multiscale::MPWideIntercommunicator>(options, intercomms);

  lMaster->RunSimulation();

  delete sharedValueBuffer;
  delete lbOrchestration;

  return (0);
}
