#include "SimulationMaster.h"

#include "usage.h"

#include <stdlib.h>

SimulationMaster::SimulationMaster(int iArgCount, char *iArgList[])
{
  mLbm = new LBM();
  mNet = new Net(iArgCount, iArgList);

  int lLbmPeriod;
  int lLbmVoxelSize;
  int lLbmSteeringSessionId;

  if (iArgCount == 8)
  {
    mMaxCycleCount = atoi(iArgList[2]);
    lLbmPeriod = atoi(iArgList[3]);
    lLbmVoxelSize = atof(iArgList[4]);
    mSnapshotsPerCycle = atoi(iArgList[5]);
    mImagesPerCycle = atoi(iArgList[6]);
    lLbmSteeringSessionId = atoi(iArgList[7]);
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
