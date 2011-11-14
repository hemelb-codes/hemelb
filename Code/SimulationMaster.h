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
#include "configuration/CommandLine.h"
#include "reporting/FileManager.h"
#include "reporting/Reporter.h"

class SimulationMaster
{
  public:
    SimulationMaster(hemelb::configuration::CommandLine &options);
    ~SimulationMaster();

    void Abort();

    bool IsCurrentProcTheIOProc();

    int GetProcessorCount();

    void RunSimulation();



  private:
    void PrintTimingData();
    void Initialise();
    void SetupReporting(); // set up the reporting file
    void PostSimulation(int iTotalTimeSteps, double iSimulationTime, bool iIsUnstable);
    unsigned int OutputPeriod(unsigned int frequency);
    void HandleActors();
    void ResetUnstableSimulation();
    void WriteLocalImages();
    void GenerateNetworkImages();
    hemelb::configuration::SimConfig *simConfig;
    hemelb::geometry::LatticeData* mLatDat;
    hemelb::reporting::FileManager* fileManager;
    typedef std::multimap<unsigned long, unsigned long> mapType;

    mapType snapshotsCompleted;
    mapType networkImagesCompleted;

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

    std::vector<hemelb::net::IteratedAction*> actors;

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
