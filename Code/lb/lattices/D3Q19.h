#ifndef HEMELB_LB_LATTICES_D3Q19_H
#define HEMELB_LB_LATTICES_D3Q19_H

#include "lb/lattices/Lattice.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      class D3Q19 : Lattice<D3Q19>
      {
        public:
          // The number of discrete velocity vectors
          static const unsigned int NUMVECTORS = 19;

          // The x, y and z components of each of the discrete velocity vectors
          static const int CX[NUMVECTORS];
          static const int CY[NUMVECTORS];
          static const int CZ[NUMVECTORS];
          static const int* discreteVelocityVectors[3];

          static const double EQMWEIGHTS[NUMVECTORS];

          // The index of the inverse direction of each discrete velocity vector
          static const int INVERSEDIRECTIONS[NUMVECTORS];
      };
    }
  }
}

#endif
