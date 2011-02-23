#ifndef HEMELB_STEERING_STEERINGCOMPONENT_H
#define HEMELB_STEERING_STEERINGCOMPONENT_H

#include "net/PhasedBroadcast.h"
#include "lb/lb.h"
#include "lb/SimulationState.h"
#include "vis/Control.h"
#include "steering/ClientConnection.h"

namespace hemelb
{
  namespace steering
  {
    class SteeringComponent : public net::PhasedBroadcast
    {
      public:
        SteeringComponent(bool* iIsNetworkSending,
                          int imagesPeriod,
                          sem_t* varLock,
                          ClientConnection* iClientConnection,
                          vis::Control* iVisControl,
                          lb::LBM* iLbm,
                          net::Net * iNet,
                          const topology::NetworkTopology *iNetTop,
                          lb::SimulationState * iSimState);

        static bool RequiresSeparateSteeringCore();

        /*
         * This function initialises all of the steering parameters, on all nodes.
         */
        void Reset();

      protected:
        void ProgressFromChildren();
        void ProgressFromParent();
        void ProgressToChildren();
        void ProgressToParent();

        void PostReceiveFromChildren();
        void PostReceiveFromParent();

        void TopNodeAction();
        void Effect();

      private:
        void AssignValues();

        const static int STEERABLE_PARAMETERS = 20;
        const static int SPREADFACTOR = 10;

        int imagesPeriod;
        bool* isNetworkSending;
        bool isConnected;
        bool updatedMouseCoords;

        sem_t* mVariableLock;
        ClientConnection* mClientConnection;
        lb::LBM* mLbm;
        lb::SimulationState* mSimState;
        vis::Control* mVisControl;
        float privateSteeringParams[STEERABLE_PARAMETERS + 1];
    };
  }
}

#endif /* HEMELB_STEERING_STEERINGCOMPONENT_H */
