#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb/lb.h"
#include "lb/StabilityTester.h"
#include "net/net.h"
#include "steering/ImageSendComponent.h"
#include "steering/SteeringComponent.h"
#include "lb/EntropyTester.h"
#include "lb/boundaries/BoundaryValues.h"
#include "util/UnitConverter.h"

class SimulationMaster
{
  public:
    SimulationMaster(int iArgCount, char *iArgList[]);
    ~SimulationMaster();

    void Abort();

    bool IsCurrentProcTheIOProc();

    int GetProcessorCount();

    void RunSimulation();



  private:
    void Initialise();
    void SetupReporting(); // set up the reporting file
    void ParseArguments(int argc, char **argv);
    void PrintUsage();
    void PostSimulation(int iTotalTimeSteps, double iSimulationTime, bool iIsUnstable);

    void PrintTimingData();

    FILE *mTimingsFile;
    std::string outputDir;
    std::string inputFile;
    std::string snapshotDirectory;
    std::string imageDirectory;

    hemelb::SimConfig *simConfig;
    hemelb::geometry::LatticeData* mLatDat;

    hemelb::steering::Network* network;
    hemelb::steering::ImageSendComponent *imageSendCpt;
    hemelb::steering::SteeringComponent* steeringCpt;

    hemelb::lb::SimulationState* mSimulationState;
    hemelb::lb::StabilityTester* mStabilityTester;
    hemelb::lb::EntropyTester* mEntropyTester;
    hemelb::lb::LBM* mLbm;
    hemelb::lb::boundaries::BoundaryValues* mInletValues;
    hemelb::lb::boundaries::BoundaryValues* mOutletValues;
    hemelb::net::Net mNet;

    hemelb::util::UnitConverter* mUnits;

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

    unsigned int snapshotsPerCycle;
    unsigned int imagesPerCycle;
    int steeringSessionId;
};

#endif /* HEMELB_SIMULATIONMASTER_H */
