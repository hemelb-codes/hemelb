
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETCOSINE_CUH
#define HEMELB_LB_IOLETS_INOUTLETCOSINE_CUH

#include "units.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetCosineGPU
      {
      public:
        distribn_t minimumSimulationDensity;
        double3 normal;
        double densityMean;
        double densityAmp;
        double phase;
        double period;
        unsigned int warmUpLength;

      public:
        __device__ distribn_t GetDensity(unsigned long timeStep) const;
      };
    }
  }
}

#endif /* HEMELB_LB_IOLETS_INOUTLETCOSINE_CUH */
