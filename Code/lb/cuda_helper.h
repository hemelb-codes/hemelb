
#ifndef HEMELB_LB_CUDA_HELPER
#define HEMELB_LB_CUDA_HELPER

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include "lb/iolets/InOutLetCosine.cuh"
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



namespace hemelb {
namespace lb {

extern iolets::InOutLetCosineGPU* inlets_dev;
extern iolets::InOutLetCosineGPU* outlets_dev;

}
}



#endif
