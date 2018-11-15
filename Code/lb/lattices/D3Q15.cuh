
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LATTICES_D3Q15_CUH
#define HEMELB_LB_LATTICES_D3Q15_CUH

#include "units.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      namespace D3Q15_GPU
      {
        __constant__ const int NUMVECTORS = 15;

        __constant__ const distribn_t CXD[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
        __constant__ const distribn_t CYD[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 };
        __constant__ const distribn_t CZD[] = { 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1 };

        __constant__ const distribn_t EQMWEIGHTS[] = {
          2.0 / 9.0,
          1.0 / 9.0,
          1.0 / 9.0,
          1.0 / 9.0,
          1.0 / 9.0,
          1.0 / 9.0,
          1.0 / 9.0,
          1.0 / 72.0,
          1.0 / 72.0,
          1.0 / 72.0,
          1.0 / 72.0,
          1.0 / 72.0,
          1.0 / 72.0,
          1.0 / 72.0,
          1.0 / 72.0
        };

        __constant__ const int INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13 };
      }
    }
  }
}

#endif /* HEMELB_LB_LATTICES_D3Q15_CUH */
