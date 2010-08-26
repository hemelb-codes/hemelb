#ifndef HEMELB_STEERING_STEERING_H
#define HEMELB_STEERING_STEERING_H

#include "steering/common/common.h"

#ifdef NO_STEER
#include "steering/off/steering.h"

#else // NO_STEER is undefined
#include "steering/on/steering.h"

#endif // NO_STEER

#endif // HEMELB_STEERING_STEERING_H
