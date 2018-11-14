
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/InOutLetCosine.cuh"



using namespace hemelb;
using namespace hemelb::lb::iolets;



__device__ distribn_t InOutLetCosineGPU::GetDensity(unsigned long timeStep) const
{
  distribn_t w = 2.0 * M_PI / period;
  distribn_t target = densityMean + densityAmp * cos(w * timeStep + phase);

  if (timeStep >= warmUpLength)
  {
    return target;
  }

  double interpolationFactor = ((double) timeStep) / ((double) warmUpLength);

  return interpolationFactor * target + (1. - interpolationFactor) * minimumSimulationDensity;
}
