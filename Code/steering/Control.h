#ifndef HEMELB_STEERING_BASIC_CONTROL_H
#define HEMELB_STEERING_BASIC_CONTROL_H

#include <pthread.h>
#include <semaphore.h>

#include "lb/SimulationState.h"
#include "lb/lb.h"
#include "vis/Control.h"
#include "steering/SteeringComponent.h"
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
                                ClientConnection* iClientConnection,
                                const lb::LbmParameters *iLbmParams);

        bool ShouldRenderForNetwork();

        sem_t nrl;

        bool is_frame_ready;
        bool sending_frame;
        Lockable<bool> isConnected;

        int send_array_length;

        bool updated_mouse_coords;

      protected:
        ~Control();

        // Is this MPI task the IO task?
        bool mIsCurrentProcTheSteeringProc;

        NetworkThread* mNetworkThread;
    };

  }

}

#endif /* HEMELB_STEERING_BASIC_CONTROL_H */
