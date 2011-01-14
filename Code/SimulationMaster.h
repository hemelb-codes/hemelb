#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb/lb.h"
#include "net.h"
#include "steering/basic/Control.h"

class SimulationMaster
{
  public:
    SimulationMaster(int iArgCount, char *iArgList[]);
    ~SimulationMaster();

    void Abort();

    bool IsCurrentProcTheIOProc();

    int GetProcessorCount();

    void RunSimulation(hemelb::SimConfig *& lSimulationConfig,
                       double iStartTime,
                       std::string image_directory,
                       std::string snapshot_directory,
                       unsigned int lSnapshotsPerCycle,
                       unsigned int lImagesPerCycle);

    void Initialise(hemelb::SimConfig *iSimConfig,
                    int iSteeringSessionid,
                    FILE * bTimingsFile);

  private:
    void PostSimulation(int iTotalTimeSteps,
                        double iSimulationTime,
                        bool iIsUnstable,
                        double iStartTime);

    void PrintTimingData();

    FILE *mTimingsFile;

    hemelb::lb::GlobalLatticeData mGlobLatDat;
    hemelb::lb::LocalLatticeData* mLocalLatDat;

    hemelb::topology::NetworkTopology* mNetworkTopology;
    hemelb::topology::TopologyManager* mTopologyManger;

    hemelb::steering::Control *steeringController;
    hemelb::lb::SimulationState mSimulationState;
    hemelb::lb::LBM *mLbm;
    Net *mNet;

    int mImagesWritten;
    int mSnapshotsWritten;

    double mDomainDecompTime;

    double mLbTime;
    double mMPISendTime;
    double mMPIWaitTime;
    double mImagingTime;
    double mSnapshotTime;
};

#endif /* HEMELB_SIMULATIONMASTER_H */
