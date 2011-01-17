#include "steering/none/Control.h"
#include "steering/common/common.h"

namespace hemelb
{

  namespace steering
  {

    Control::Control(bool isCurrentProcTheSteeringProc)
    {
    }

    Control::~Control()
    {
    }

    // Kick off the networking thread
    void Control::StartNetworkThread(lb::LBM* lbm,
                                     lb::SimulationState *iSimState,
                                     const lb::LbmParameters *iLbmParams)
    {
    }

    // Broadcast the steerable parameters to all tasks.
    void Control::UpdateSteerableParameters(bool shouldRenderForSnapshot,
                                            int* perform_rendering,
                                            hemelb::lb::SimulationState &iSimulationState,
                                            hemelb::vis::Control* visController,
                                            lb::LBM* lbm)
    {
    }

    // Do we need to render a frame for the client?
    bool Control::ShouldRenderForNetwork()
    {
      return false;
    }

  }

}
