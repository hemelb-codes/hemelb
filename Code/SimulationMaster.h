#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb/lb.hpp"
#include "lb/StabilityTester.h"
#include "net/net.h"
#include "steering/ImageSendComponent.h"
#include "steering/SteeringComponent.h"
#include "lb/EntropyTester.h"
#include "lb/boundaries/BoundaryValues.h"
#include "util/UnitConverter.h"
#include "configuration/CommandLine.h"
#include "io/PathManager.h"
#include "reporting/Reporter.h"
#include "reporting/Timers.h"
#include "lb/IncompressibilityChecker.hpp"

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
    typedef hemelb::D3Q15 latticeType;

    void Initialise();
    void SetupReporting(); // set up the reporting file
    unsigned int OutputPeriod(unsigned int frequency);
    void HandleActors();
    void ResetUnstableSimulation();
    void WriteLocalImages();
    void GenerateNetworkImages();
    hemelb::configuration::SimConfig *simConfig;
    hemelb::geometry::LatticeData* latticeData;
    hemelb::io::PathManager* fileManager;
    hemelb::reporting::Timers timings;
    hemelb::reporting::Reporter* reporter;
    typedef std::multimap<unsigned long, unsigned long> MapType;

    MapType snapshotsCompleted;
    MapType networkImagesCompleted;

    hemelb::steering::Network* network;
    hemelb::steering::ImageSendComponent *imageSendCpt;
    hemelb::steering::SteeringComponent* steeringCpt;

    hemelb::lb::SimulationState* simulationState;
    hemelb::lb::StabilityTester* stabilityTester;
    hemelb::lb::EntropyTester<latticeType>* entropyTester;

    /** Actor in charge of checking the maximum density difference across the domain */
    hemelb::lb::IncompressibilityChecker<>* incompressibilityChecker;

    hemelb::lb::LBM<latticeType>* latticeBoltzmannModel;
    hemelb::lb::boundaries::BoundaryValues* inletValues;
    hemelb::lb::boundaries::BoundaryValues* outletValues;
    hemelb::net::Net communicationNet;

    hemelb::util::UnitConverter* unitConvertor;

    hemelb::vis::Control* visualisationControl;

    std::vector<hemelb::net::IteratedAction*> actors;

    unsigned int snapshotsPerCycle;
    unsigned int imagesPerCycle;
    int steeringSessionId;
};

#endif /* HEMELB_SIMULATIONMASTER_H */
