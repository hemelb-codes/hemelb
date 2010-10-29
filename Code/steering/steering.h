#ifndef HEMELB_STEERING_STEERING_H
#define HEMELB_STEERING_STEERING_H

#include "steering/common/common.h"

#ifdef HEMELB_STEERING_LIB_basic
  #include "steering/basic/Control.h"

#else
  #include "steering/none/steering.h"

#endif

#endif // HEMELB_STEERING_STEERING_H
