
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/lattices/D3Q7.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<>
      lb::lattices::LatticeInfo* lb::lattices::Lattice<D3Q7>::singletonInfo = NULL;
           
      const int D3Q7::CX[] = { 0, 1, -1, 0, 0, 0, 0 };
      const int D3Q7::CY[] = { 0, 0, 0, 1, -1, 0, 0 };
      const int D3Q7::CZ[] = { 0, 0, 0, 0, 0, 1, -1 };
      const int* D3Q7::discreteVelocityVectors[] = { CX, CY, CZ };

      const distribn_t D3Q7::CXD[] = { 0.0, 1.0, -1.0, 0.0,  0.0, 0.0,  0.0};
      const distribn_t D3Q7::CYD[] = { 0.0, 0.0,  0.0, 1.0, -1.0, 0.0,  0.0};
      const distribn_t D3Q7::CZD[] = { 0.0, 0.0,  0.0, 0.0,  0.0, 1.0, -1.0};            
      
      const distribn_t D3Q7::EQMWEIGHTS[] = { 1.0 / 4.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0 };

      const Direction D3Q7::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5 };
    }
  }
}
