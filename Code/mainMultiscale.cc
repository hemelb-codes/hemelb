#include "configuration/CommandLine.h"
#include "SimulationMaster.h"
#include <multiscale/MultiscaleSimulationMaster.h>
#include <resources/Resource.h>
#include <multiscale/Intercommunicator.h>
#include "multiscale/mpwide/MPWideIntercommunicator.h"

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

  std::map<std::string,double> *pbuffer;
  std::map<std::string,bool> *orchestrationLB;

  hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc,argv);
  
  hemelb::multiscale::MultiscaleSimulationMaster<hemelb::multiscale::MPWideIntercommunicator> *lMaster;

  pbuffer = new std::map<std::string, double>();
  std::map<std::string, double> &buffer = *pbuffer;

  orchestrationLB=new std::map<std::string,bool>();
  std::map<std::string,bool> &rorchestrationLB=*orchestrationLB;
  rorchestrationLB["boundary1_pressure"] = false;
  rorchestrationLB["boundary2_pressure"] = false;
  rorchestrationLB["boundary1_velocity"] = true;
  rorchestrationLB["boundary2_velocity"] = true;

  hemelb::multiscale::mpwide::mpwide_config_file = "../../config_files/MPWSettings.cfg";
  hemelb::multiscale::MPWideIntercommunicator intercomms(*pbuffer,*orchestrationLB);
  //TODO: Add an IntercommunicatorImplementation?

  lMaster = new hemelb::multiscale::MultiscaleSimulationMaster<hemelb::multiscale::MPWideIntercommunicator>(options, intercomms);

  lMaster->RunSimulation();

  return (0);
}
