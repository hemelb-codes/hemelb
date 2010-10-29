#ifndef HEMELB_STEERING_BASIC_STEERING_H
#define HEMELB_STEERING_BASIC_STEERING_H

#include <pthread.h>
#include <semaphore.h>

#include "lb.h"
#include "vis/Control.h"
#include "steering/basic/Control.h"

namespace hemelb
{
  namespace steering
  {

    void  *hemeLB_steer (void*);
    void* hemeLB_network(void *ptr);
  }
}


#endif // HEME_STEERING_BASIC_STEERING_H
