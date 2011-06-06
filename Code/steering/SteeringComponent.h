#ifndef HEMELB_STEERING_STEERINGCOMPONENT_H
#define HEMELB_STEERING_STEERINGCOMPONENT_H

#include "net/PhasedBroadcastRegular.h"
#include "lb/lb.h"
#include "lb/SimulationState.h"
#include "vis/DomainStats.h"
#include "vis/Control.h"
#include "steering/Network.h"

namespace hemelb
{
  namespace steering
  {
    class SteeringComponent : public net::PhasedBroadcastRegular<false, 1, 0, true, false>
    {
      public:
        SteeringComponent(Network* iNetwork,
                          vis::Control* iVisControl,
                          lb::LBM* iLbm,
                          net::Net * iNet,
                          lb::SimulationState * iSimState);

        static bool RequiresSeparateSteeringCore();

        /*
         * This function initialises all of the steering parameters, on all nodes.
         */
        void Reset();

        bool readyForNextImage;
        bool updatedMouseCoords;

      protected:
        void ProgressFromParent(unsigned long splayNumber);
        void ProgressToChildren(unsigned long splayNumber);

        void TopNodeAction();
        void Effect();

      private:
        void AssignValues();

        const static int STEERABLE_PARAMETERS = 20;
        const static unsigned int SPREADFACTOR = 10;

        bool isConnected;

        Network* mNetwork;
        lb::LBM* mLbm;
        lb::SimulationState* mSimState;
        vis::Control* mVisControl;
        float privateSteeringParams[STEERABLE_PARAMETERS + 1];
    };
  }
}

#endif /* HEMELB_STEERING_STEERINGCOMPONENT_H */
