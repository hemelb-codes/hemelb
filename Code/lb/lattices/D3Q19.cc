
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/lattices/D3Q19.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      const int D3Q19::CX[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0 };
      const int D3Q19::CY[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1 };
      const int D3Q19::CZ[] = { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1 };

      const distribn_t D3Q19::CXD[] = { 0.0, 1.0, -1.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0,  1.0, 
                                        -1.0, 1.0, -1.0,  1.0, -1.0, 0.0,  0.0,  0.0,  0.0 };
      const distribn_t D3Q19::CYD[] = { 0.0, 0.0,  0.0, 1.0, -1.0, 0.0,  0.0, 1.0, -1.0, -1.0,
                                        1.0, 0.0,  0.0,  0.0,  0.0, 1.0, -1.0,  1.0, -1.0 };
      const distribn_t D3Q19::CZD[] = { 0.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0, 0.0,  0.0,  0.0,
                                        0.0, 1.0, -1.0, -1.0,  1.0, 1.0, -1.0, -1.0,  1.0 };

      const int* D3Q19::discreteVelocityVectors[] = { CX, CY, CZ };

      const distribn_t D3Q19::EQMWEIGHTS[] = { 1.0 / 3.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
                                               1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
                                               1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
                                               1.0 / 36.0 };
                  
      const Direction D3Q19::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };
    }
  }
}
