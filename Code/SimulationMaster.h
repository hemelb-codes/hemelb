#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb/lb.h"
#include "lb/StabilityTester.h"
#include "net/net.h"
#include "steering/ImageSendComponent.h"
#include "steering/SteeringComponent.h"
#include "lb/EntropyTester.h"

class SimulationMaster
{
  public:
    SimulationMaster(int iArgCount, char *iArgList[]);
    ~SimulationMaster();

    void Abort();

    bool IsCurrentProcTheIOProc();

    int GetProcessorCount();

    void RunSimulation(std::string image_directory,
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
                        bool iIsUnstable);

    void PrintTimingData();

    FILE *mTimingsFile;

    hemelb::geometry::LatticeData* mLatDat;

    hemelb::steering::Network* network;
    hemelb::steering::ImageSendComponent *imageSendCpt;
    hemelb::steering::SteeringComponent* steeringCpt;

    hemelb::lb::SimulationState* mSimulationState;
    hemelb::lb::StabilityTester* mStabilityTester;
    hemelb::lb::EntropyTester* mEntropyTester;
    hemelb::lb::LBM *mLbm;
    hemelb::net::Net mNet;

    hemelb::vis::Control* mVisControl;

    int mImagesWritten;
    int mSnapshotsWritten;

    double mCreationTime;

    double mSnapshotTime;
    double mDomainDecompTime;
    double mFileReadTime;
    double mNetInitialiseTime;

    double mMPISendTime;
    double mMPIWaitTime;
};

#endif /* HEMELB_SIMULATIONMASTER_H */
