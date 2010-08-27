#ifndef HEMELB_STEERING_STEERING_H
#define HEMELB_STEERING_STEERING_H

#include "steering/common/common.h"

#if HEMELB_STEERING_LIB == basic
  #include "steering/basic/steering.h"

#else
  #include "steering/none/steering.h"

#endif

#endif // HEMELB_STEERING_STEERING_H
