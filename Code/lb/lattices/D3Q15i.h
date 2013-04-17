//
// Copyright (C) University College London, 2007-2013, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

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

          static const double EQMWEIGHTS[NUMVECTORS];

          // The index of the inverse direction of each discrete velocity vector
          static const Direction INVERSEDIRECTIONS[NUMVECTORS];
      };

    } /* namespace lattices */
  } /* namespace lb */
} /* namespace hemelb */
#endif /* HEMELB_LB_LATTICES_D3Q15I_H */
