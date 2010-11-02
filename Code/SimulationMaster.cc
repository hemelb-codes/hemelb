#include "SimulationMaster.h"
#include "SimConfig.h"

#include "usage.h"
#include "debug/Debugger.h"

#include <stdlib.h>

SimulationMaster::SimulationMaster(int iArgCount, char *iArgList[])
{
  //TODO SLOW INITIALISATION. Make these statically allocated. Currently no reason not to.
  mLbm = new LBM();
  mNet = new Net(iArgCount, iArgList);

  hemelb::debug::Debugger::Init(iArgList[0]);
}

SimulationMaster::~SimulationMaster()
{
  delete mNet;
  delete mLbm;
}

LBM *SimulationMaster::GetLBM()
{
  return mLbm;
}

Net *SimulationMaster::GetNet()
{
  return mNet;
}
