#ifndef HEMELB_SIMULATIONMASTER_H
#define HEMELB_SIMULATIONMASTER_H

#include "lb.h"
#include "net.h"

class SimulationMaster
{
  public:
    SimulationMaster(int iArgCount, char *iArgList[]);
    ~SimulationMaster();

    // TODO: These functions are a temporary hack for moving stuff into this class from the main function.
    LBM *GetLBM();
    Net *GetNet();

  private:

    LBM *mLbm;
    Net *mNet;
};

#endif /* HEMELB_SIMULATIONMASTER_H */
