// In this file all the declarations and some simple functions are
// reported.
// The structs are defined here.

// Global coordinate means coordinate within the entire system, not the
// coordinate on one proc.
#ifndef __config_h__
#define __config_h__

#include "mpiInclude.h"

// These 3 lines break things on many machines
// #ifndef int64_t
// #define int64_t long int
// #endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "constants.h"
#include "lb.h"
#include "rt.h"

#ifndef NO_STEER
#include <pthread.h>
#include <semaphore.h>
#endif





#endif                  // __config_h__
