
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/lattices/D3Q27ad.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<>
      lb::lattices::LatticeInfo* lb::lattices::Lattice<D3Q27ad>::singletonInfo = NULL;

      const int D3Q27ad::CX[] =
          { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
      const int D3Q27ad::CY[] =
          { 0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1 };
      const int D3Q27ad::CZ[] =
          { 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1 };
      
      const distribn_t  D3Q27ad::CXD[] = 
          { 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0,  -1.0, 1.0, -1.0, 0.0, 0.0, 
            0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };
      const distribn_t  D3Q27ad::CYD[] = 
          { 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0,
           -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0 };
      const distribn_t  D3Q27ad::CZD[] = 
         { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,  -1.0, 1.0, 1.0, -1.0, -1.0,
          1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0 };
      
      const int* D3Q27ad::discreteVelocityVectors[] = { CX, CY, CZ };

      const distribn_t D3Q27ad::EQMWEIGHTS[] = { 8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,
                                               2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
                                               1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
                                               1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0
                                                   / 216.0,
                                               1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0 };

      const Direction D3Q27ad::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17,
                                                     20, 19, 22, 21, 24, 23, 26, 25 };
    }
  }
}
