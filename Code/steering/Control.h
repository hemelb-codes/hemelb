#ifndef HEMELB_STEERING_BASIC_CONTROL_H
#define HEMELB_STEERING_BASIC_CONTROL_H

#include <pthread.h>
#include <semaphore.h>

#include "lb/SimulationState.h"
#include "lb/lb.h"
#include "vis/Control.h"
#include "steering/basic/NetworkThread.h"
#include "steering/basic/Lockable.h"

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

        void UpdateSteerableParameters(bool shouldRenderForSnapshot,
                                       hemelb::lb::SimulationState &iSimulationState,
                                       hemelb::vis::Control* visController,
                                       lb::LBM* lbm);
        bool ShouldRenderForNetwork();

        bool RequiresSeparateSteeringCore() const;

        sem_t nrl;
        sem_t steering_var_lock;

        bool is_frame_ready;
        bool sending_frame;
        Lockable<bool> isConnected;

        int send_array_length;

        bool updated_mouse_coords;

      protected:
        ~Control();

        // Is this MPI task the IO task?
        bool mIsCurrentProcTheSteeringProc;

        // Do the MPI send to spread the params
        void
        BroadcastSteerableParameters(hemelb::lb::SimulationState &lSimulationState,
                                     vis::Control *visControl,
                                     lb::LBM* lbm);

        NetworkThread* mNetworkThread;
    };

  }

}

#endif /* HEMELB_STEERING_BASIC_CONTROL_H */
