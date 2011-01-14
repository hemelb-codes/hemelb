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

    Control* Control::Init(bool isCurrentProcTheSteeringProc)
    {
      /* Static member function that implements the singleton pattern.
       */
      if (Control::singleton == NULL)
      {
        Control::singleton = new Control(isCurrentProcTheSteeringProc);
      }
      return Control::singleton;
    }

    Control* Control::Get(void)
    {
      // Get the single instance.
      return Control::singleton;
    }

    // Init static members
    Control* Control::singleton = NULL;

    // Kick off the networking thread
    void Control::StartNetworkThread(LBM* lbm,
                                     lb::SimulationState *iSimState,
                                     const lb::LbmParameters *iLbmParams)
    {
    }

    // Broadcast the steerable parameters to all tasks.
    void Control::UpdateSteerableParameters(bool shouldRenderForSnapshot,
                                            int* perform_rendering,
                                            hemelb::lb::SimulationState &iSimulationState,
                                            hemelb::vis::Control* visController,
                                            LBM* lbm)
    {
    }

    // Do we need to render a frame for the client?
    bool Control::ShouldRenderForNetwork()
    {
      return false;
    }

  }

}
