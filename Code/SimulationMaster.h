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

    void RunSimulation(FILE *iTimingsFile,
                       hemelb::SimConfig *& lSimulationConfig,
                       double iStartTime,
                       char image_directory[256],
                       char snapshot_directory[256],
                       unsigned int lSnapshotsPerCycle,
                       unsigned int lImagesPerCycle);

    void Initialise(hemelb::SimConfig *iSimConfig,
                    int iSteeringSessionid,
                    FILE *bTimingsFile);

  private:
    void PostSimulation(int iTotalTimeSteps,
                        double iSimulationTime,
                        FILE *timings_ptr,
                        bool iIsUnstable,
                        double iStartTime);

    hemelb::steering::Control *steeringController;
    hemelb::lb::SimulationState mSimulationState;
    LBM *mLbm;
    Net *mNet;
};

#endif /* HEMELB_SIMULATIONMASTER_H */
