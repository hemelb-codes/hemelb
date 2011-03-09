#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb/lb.h"
#include "lb/StabilityTester.h"
#include "net/net.h"
#include "steering/ImageSendComponent.h"
#include "steering/SteeringComponent.h"

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
                    unsigned int iImagesPerCycle,
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

    hemelb::steering::ClientConnection* clientConnection;
    hemelb::steering::ImageSendComponent *imageSendCpt;
    hemelb::steering::SteeringComponent* steeringCpt;

    hemelb::lb::SimulationState mSimulationState;
    hemelb::lb::StabilityTester* mStabilityTester;
    hemelb::lb::LBM *mLbm;
    hemelb::net::Net *mNet;

    hemelb::vis::Control* mVisControl;

    int mImagesWritten;
    int mSnapshotsWritten;

    double mDomainDecompTime;
    double mFileReadTime;
    double mNetInitialiseTime;

    double mLbTime;
    double mMPISendTime;
    double mMPIWaitTime;
    double mImagingTime;
    double mSnapshotTime;
};

#endif /* HEMELB_SIMULATIONMASTER_H */
