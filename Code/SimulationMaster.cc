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

  if (iArgCount == 8)
  {
    mMaxCycleCount = atoi(iArgList[2]);
    mLbm->period = atoi(iArgList[3]);
    mLbm->voxel_size = atof(iArgList[4]);
    mSnapshotsPerCycle = atoi(iArgList[5]);
    mImagesPerCycle = atoi(iArgList[6]);
    mLbm->steering_session_id = atoi(iArgList[7]);
  }
  else
  {
    if (mNet->IsCurrentProcTheIOProc())
    {
      Usage::printUsage(iArgList[0]);
    }

    mNet->Abort();
  }
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
