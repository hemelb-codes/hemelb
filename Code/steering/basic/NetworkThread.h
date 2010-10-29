#ifndef HEMELB_STEERING_BASIC_NETWORKTHREAD
#define HEMELB_STEERING_BASIC_NETWORKTHREAD

#include "steering/basic/Threadable.h"

#include "lb.h"

namespace hemelb
{
  namespace steering
  {
    class Control;

    class NetworkThread : public Threadable
    {
      public:
        NetworkThread(LBM* lbm, Control* steeringController);
      private:
        void DoWork(void);

        LBM* mLbm;
        Control* mSteeringController;

        pthread_attr_t* GetPthreadAttributes(void);
        double frameTiming(void);

        static const unsigned int MYPORT;
        static const unsigned int CONNECTION_BACKLOG;
    };

  } // namespace steering
} //namespace hemelb

#endif // HEMELB_STEERING_BASIC_NETWORKTHREAD
