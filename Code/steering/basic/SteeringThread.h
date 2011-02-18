#ifndef HEMELB_STEERING_BASIC_STEERINGTHREAD
#define HEMELB_STEERING_BASIC_STEERINGTHREAD

#include <semaphore.h>
#include "steering/basic/ClientConnection.h"
#include "steering/basic/Threadable.h"

namespace hemelb
{
  namespace steering
  {

    class SteeringThread : public Threadable
    {
      public:
        SteeringThread(ClientConnection* iClientConn, sem_t * bVariableEditSempahore);

      private:
        void DoWork(void);

        ClientConnection* mClientConn;
        sem_t* mSteeringVariableLock;
    };

  }
}
#endif // HEMELB_STEERING_BASIC_STEERINGTHREAD
