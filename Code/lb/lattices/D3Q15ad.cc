
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/lattices/D3Q15ad.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<>
      lb::lattices::LatticeInfo* lb::lattices::Lattice<D3Q15ad>::singletonInfo = NULL;

      const int D3Q15ad::CX[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
      const int D3Q15ad::CY[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 };
      const int D3Q15ad::CZ[] = { 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1 };
      const int* D3Q15ad::discreteVelocityVectors[] = { CX, CY, CZ };
      
      const distribn_t D3Q15ad::CXD[] = { 0.0, 1.0, -1.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0 };
      const distribn_t D3Q15ad::CYD[] = { 0.0, 0.0,  0.0, 1.0, -1.0, 0.0,  0.0, 1.0, -1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0 };
      const distribn_t D3Q15ad::CZD[] = { 0.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0, 1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0 };

      const distribn_t D3Q15ad::EQMWEIGHTS[] = { 2.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0
                                                    / 9.0,
                                                1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0,
                                                1.0 / 72.0, 1.0 / 72.0 };

      const Direction D3Q15ad::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13 };
    } /* namespace lattices */
  } /* namespace lb */
} /* namespace hemelb */
