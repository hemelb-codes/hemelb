#ifndef HEMELB_STEERING_BASIC_NETWORKTHREAD
#define HEMELB_STEERING_BASIC_NETWORKTHREAD

#ifndef NO_STEER
#include <pthread.h>
#endif

#include "steering/basic/Threadable.h"
#include "lb/SimulationState.h"
#include "lb/lb.h"

namespace hemelb
{
  namespace steering
  {
    class Control;

    class NetworkThread : public Threadable
    {
      public:
        NetworkThread(lb::LBM* lbm,
                      Control* steeringController,
                      lb::SimulationState* iSimState,
                      const lb::LbmParameters* iLbmParams);

        ~NetworkThread();

#ifndef NO_STEER
        static pthread_mutex_t var_lock;
#endif

      private:
        void DoWork(void);

        lb::LBM* mLbm;
        Control* mSteeringController;
        lb::SimulationState* mSimState;
        const lb::LbmParameters* mLbmParams;

        pthread_attr_t* GetPthreadAttributes(void);
        double frameTiming(void);

        static const unsigned int MYPORT = 65250;
        static const unsigned int CONNECTION_BACKLOG = 10;

        // data per pixel inlude the data for the pixel location and 4 colours
        //in RGB format, thus #bytes per pixel are (sizeof(int)+4
        //rgb)=sizeof(int)+4*3*sizeof(unsigned char))
        static const int bytes_per_pixel_data = sizeof(int) + 4 * sizeof(unsigned char);
        // one int for colour_id and one for pixel id
        static const u_int pixel_data_bytes = COLOURED_PIXELS_MAX * bytes_per_pixel_data;
        // it is assumed that the frame size is the only detail
        static const u_int frame_details_bytes = 1 * sizeof(int);
        char* xdrSendBuffer_pixel_data;
        char* xdrSendBuffer_frame_details;

        void setRenderState(int val);
    };

  } // namespace steering
} //namespace hemelb

#endif // HEMELB_STEERING_BASIC_NETWORKTHREAD
