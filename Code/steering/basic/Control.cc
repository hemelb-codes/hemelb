/*
 * Control.cc
 *
 *  Created on: Oct 27, 2010
 *      Author: rupert
 */
#include <pthread.h>
#include <semaphore.h>

#include "lb/lb.h"
#include "steering/Control.h"
#include "steering/ClientConnection.h"

namespace hemelb
{

  namespace steering
  {

    Control::Control(bool isCurrentProcTheSteeringProc) :
      is_frame_ready(false), sending_frame(false), isConnected(false), updated_mouse_coords(false),
          mIsCurrentProcTheSteeringProc(isCurrentProcTheSteeringProc)
    {
      sem_init(&nrl, 0, 1);
    }

    Control::~Control()
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        delete mNetworkThread;
      }
      sem_destroy(&nrl);
    }

    // Kick off the networking thread
    void Control::StartNetworkThread(lb::LBM* lbm,
                                     lb::SimulationState *iSimState,
                                     ClientConnection* iClientConnection,
                                     const lb::LbmParameters *iLbmParams)
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        mNetworkThread = new NetworkThread(lbm, this, iSimState, iLbmParams, iClientConnection);
        mNetworkThread->Run();
      }
    }

    // Do we need to render a frame for the client?
    bool Control::ShouldRenderForNetwork()
    {
      // Iff we're connected and not sending a
      // frame, we'd better render a new one!
      return isConnected.GetValue() && !sending_frame;
    }
  }
}
