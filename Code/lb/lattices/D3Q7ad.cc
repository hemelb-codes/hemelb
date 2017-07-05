
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/lattices/D3Q7ad.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<>
      lb::lattices::LatticeInfo* lb::lattices::Lattice<D3Q7ad>::singletonInfo = NULL;

      const int D3Q7ad::CX[] = { 0, 1, -1, 0, 0, 0, 0 };
      const int D3Q7ad::CY[] = { 0, 0, 0, 1, -1, 0, 0 };
      const int D3Q7ad::CZ[] = { 0, 0, 0, 0, 0, 1, -1 };
      const int* D3Q7ad::discreteVelocityVectors[] = { CX, CY, CZ };
      
      const distribn_t D3Q7ad::CXD[] = { 0.0, 1.0, -1.0, 0.0,  0.0, 0.0,  0.0 };
      const distribn_t D3Q7ad::CYD[] = { 0.0, 0.0,  0.0, 1.0, -1.0, 0.0,  0.0 };
      const distribn_t D3Q7ad::CZD[] = { 0.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0 };

      const distribn_t D3Q7ad::EQMWEIGHTS[] = { 1.0 / 4.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0 };

      const Direction D3Q7ad::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5 };
    } /* namespace lattices */
  } /* namespace lb */
} /* namespace hemelb */
