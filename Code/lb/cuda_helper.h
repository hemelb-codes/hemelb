
#ifndef HEMELB_LB_CUDA_HELPER
#define HEMELB_LB_CUDA_HELPER

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include "units.h"



#define CUDA_SAFE_CALL(x)                         \
{                                                 \
  cudaError_t error = x;                          \
  if ( error != cudaSuccess ) {                   \
    const char *name = cudaGetErrorName(error);   \
    const char *str = cudaGetErrorString(error);  \
    std::cerr << "\n"                             \
              << "CUDA Error at " #x "\n"         \
              << name << ": " << str << "\n";     \
    exit(1);                                      \
  }                                               \
}



typedef struct
{
  double minimumSimulationDensity;
  double3 normal;
  double densityMean;
  double densityAmp;
  double phase;
  double period;
  unsigned int warmUpLength;
} iolet_cosine_t;



extern iolet_cosine_t* inlets_dev;
extern iolet_cosine_t* outlets_dev;



#endif
