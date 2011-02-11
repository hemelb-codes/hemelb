#ifndef HEMELB_STEERING_NONE_CONTROL_H
#define HEMELB_STEERING_NONE_CONTROL_H

#include "lb/SimulationState.h"
#include "lb/lb.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace steering
  {

    class Control
    {
      public:
        Control(bool isCurrentProcTheSteeringProc);

        void StartNetworkThread(lb::LBM* lbm,
                                lb::SimulationState *iSimState,
                                const lb::LbmParameters *iLbmParams);

        void
        UpdateSteerableParameters(bool shouldRenderForSnapshot,
                                  hemelb::lb::SimulationState &iSimulationState,
                                  hemelb::vis::Control* visController,
                                  lb::LBM* lbm);
        bool ShouldRenderForNetwork();

        bool RequiresSeparateSteeringCore() const;

      protected:
        ~Control();
    };

  }

}

#endif /* HEMELB_STEERING_NONE_CONTROL_H */
