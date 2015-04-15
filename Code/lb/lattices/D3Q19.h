// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_LATTICES_D3Q19_H
#define HEMELB_LB_LATTICES_D3Q19_H

#include "lb/lattices/Lattice.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      class D3Q19 : public Lattice<D3Q19>
      {
        public:
          // The number of discrete velocity vectors
          static const Direction NUMVECTORS = 19;

          // The x, y and z components of each of the discrete velocity vectors
          static const int CX[NUMVECTORS];
          static const int CY[NUMVECTORS];
          static const int CZ[NUMVECTORS];

          // the same in double (in order to prevent int->double conversions), and aligned to 16B
          static const distribn_t CXD[NUMVECTORS] __attribute__((aligned(16)));
          static const distribn_t CYD[NUMVECTORS] __attribute__((aligned(16)));
          static const distribn_t CZD[NUMVECTORS] __attribute__((aligned(16)));

          static const int* discreteVelocityVectors[3];

          static const double EQMWEIGHTS[NUMVECTORS] __attribute__((aligned(16)));

          // The index of the inverse direction of each discrete velocity vector
          static const Direction INVERSEDIRECTIONS[NUMVECTORS];
      };
    }
  }
}

#endif
