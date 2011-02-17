#ifndef HEMELB_STEERING_BASIC_STEERINGTHREAD
#define HEMELB_STEERING_BASIC_STEERINGTHREAD

#include <semaphore.h>
#include "steering/basic/Threadable.h"

namespace hemelb
{
  namespace steering
  {

    class SteeringThread : public Threadable
    {
      public:
        SteeringThread(int fd, sem_t * bVariableEditSempahore);

      private:
        void DoWork(void);

        sem_t* mSteeringVariableLock;
        int mSocketFileDescriptor;
    };

  }
}
#endif // HEMELB_STEERING_BASIC_STEERINGTHREAD
