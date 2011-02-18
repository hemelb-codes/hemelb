#ifndef HEMELB_STEERING_BASIC_NETWORKTHREAD
#define HEMELB_STEERING_BASIC_NETWORKTHREAD

#include <pthread.h>

#include "steering/basic/ClientConnection.h"
#include "steering/basic/Threadable.h"
#include "lb/SimulationState.h"
#include "lb/lb.h"

namespace hemelb
{
  namespace steering
  {
    class Control;

    /**
     * Class to handle the business of sending imaging and steering data over a network connection,
     * run on a separate thread.
     */
    class NetworkThread : public Threadable
    {
      public:
        NetworkThread(lb::LBM* lbm,
                      Control* steeringController,
                      lb::SimulationState* iSimState,
                      const lb::LbmParameters* iLbmParams,
                      ClientConnection* iClientConnection);

        ~NetworkThread();

      private:
        ClientConnection* mClientConnection;

        static pthread_mutex_t var_lock;

        void DoWork(void);

        void HandleBrokenPipe();

        lb::LBM* mLbm;
        Control* mSteeringController;
        lb::SimulationState* mSimState;
        const lb::LbmParameters* mLbmParams;

        pthread_attr_t* GetPthreadAttributes(void);
        double frameTiming(void);

        char* xdrSendBuffer_pixel_data;

        // data per pixel inlude the data for the pixel location and 4 colours
        //in RGB format, thus #bytes per pixel are (sizeof(int)+4
        //rgb)=sizeof(int)+4*3*sizeof(unsigned char))
        static const int bytes_per_pixel_data = sizeof(int) + 4 * sizeof(unsigned char);
        // one int for colour_id and one for pixel id
        static const u_int pixel_data_bytes = COLOURED_PIXELS_MAX * bytes_per_pixel_data;
        // it is assumed that the frame size is the only detail
        static const u_int frame_details_bytes = 1 * sizeof(int);

        void setRenderState(int val);
    };

  } // namespace steering
} //namespace hemelb

#endif // HEMELB_STEERING_BASIC_NETWORKTHREAD
