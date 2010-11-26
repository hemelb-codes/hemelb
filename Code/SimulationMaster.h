#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb.h"
#include "net.h"
#include "steering/basic/Control.h"

class SimulationMaster
{
  public:
    SimulationMaster(int iArgCount, char *iArgList[]);
    ~SimulationMaster();

    // TODO: These functions are a temporary hack for moving stuff into this class from the main function.
    LBM *GetLBM();
    Net *GetNet();

    void RunSimulation(hemelb::SimConfig *& lSimulationConfig,
                       double iStartTime,
                       std::string image_directory,
                       std::string snapshot_directory,
                       unsigned int lSnapshotsPerCycle,
                       unsigned int lImagesPerCycle);

    void Initialise(hemelb::SimConfig *iSimConfig,
                    int iSteeringSessionid,
                    FILE *bTimingsFile);

    // TODO Temporary hack while refactoring, so the main method can find out if it's
    // on the IO rank or not.
    hemelb::topology::NetworkTopology mNetworkTopology;

  private:
    void PostSimulation(int iTotalTimeSteps,
                        double iSimulationTime,
                        bool iIsUnstable,
                        double iStartTime);

    void PrintTimingData(int iSignal);

    FILE *mTimingsFile;

    hemelb::lb::GlobalLatticeData mGlobLatDat;
    hemelb::lb::LocalLatticeData mLocalLatDat;

    hemelb::topology::TopologyManager mTopologyManger;

    hemelb::steering::Control *steeringController;
    hemelb::lb::SimulationState mSimulationState;
    LBM *mLbm;
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
