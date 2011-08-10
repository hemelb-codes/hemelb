#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb/lb.h"
#include "lb/StabilityTester.h"
#include "net/net.h"
#include "steering/ImageSendComponent.h"
#include "steering/SteeringComponent.h"
#include "lb/EntropyTester.h"
#include "lb/BoundaryComms.h"
#include "util/UnitConverter.h"

class SimulationMaster
{
  public:
    SimulationMaster(int iArgCount, char *iArgList[]);
    ~SimulationMaster();

    void Abort();

    bool IsCurrentProcTheIOProc();

    int GetProcessorCount();

    void RunSimulation(double iStartTime,
                       std::string image_directory,
                       std::string snapshot_directory,
                       unsigned int lSnapshotsPerCycle,
                       unsigned int lImagesPerCycle,
                       bool doEntropyTest);

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

    hemelb::geometry::LatticeData* mLatDat;

    hemelb::steering::Network* network;
    hemelb::steering::ImageSendComponent *imageSendCpt;
    hemelb::steering::SteeringComponent* steeringCpt;

    hemelb::lb::SimulationState* mSimulationState;
    hemelb::lb::StabilityTester* mStabilityTester;
    hemelb::lb::EntropyTester* mEntropyTester;
    hemelb::lb::LBM* mLbm;
    hemelb::lb::BoundaryComms* mBoundaryComms;
    hemelb::net::Net mNet;

    hemelb::util::UnitConverter* mUnits;

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
