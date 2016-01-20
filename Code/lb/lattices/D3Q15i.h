
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_LATTICES_D3Q15I_H
#define HEMELB_LB_LATTICES_D3Q15I_H

#include "lb/lattices/IncompressibleLattice.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {

      class D3Q15i : public hemelb::lb::lattices::IncompressibleLattice<D3Q15i>
      {
        public:
          // The number of discrete velocity vectors
          static const Direction NUMVECTORS = 15;

          // The x, y and z components of each of the discrete velocity vectors
          static const int CX[NUMVECTORS];
          static const int CY[NUMVECTORS];
          static const int CZ[NUMVECTORS];
          static const int* discreteVelocityVectors[3];
          
          // the same in double (in order to prevent int->double conversions), and aligned to 16B
          static const distribn_t CXD[NUMVECTORS] __attribute__((aligned(16)));
          static const distribn_t CYD[NUMVECTORS] __attribute__((aligned(16)));
          static const distribn_t CZD[NUMVECTORS] __attribute__((aligned(16)));

          static const double EQMWEIGHTS[NUMVECTORS] __attribute__((aligned(16)));

          // The index of the inverse direction of each discrete velocity vector
          static const Direction INVERSEDIRECTIONS[NUMVECTORS] ;
      };

    } /* namespace lattices */
  } /* namespace lb */
} /* namespace hemelb */
#endif /* HEMELB_LB_LATTICES_D3Q15I_H */
