#ifndef HEMELB_LB_LATTICES_D3Q27_H
#define HEMELB_LB_LATTICES_D3Q27_H

#include "lb/lattices/Lattice.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      class D3Q27 : public Lattice<D3Q27>
      {
        public:
          // The number of discrete velocity vectors
          static const Direction NUMVECTORS = 27;

          // The x, y and z components of each of the discrete velocity vectors
          static const int CX[NUMVECTORS];
          static const int CY[NUMVECTORS];
          static const int CZ[NUMVECTORS];
          static const int* discreteVelocityVectors[3];

          static const double EQMWEIGHTS[NUMVECTORS];

          // The index of the inverse direction of each discrete velocity vector
          static const Direction INVERSEDIRECTIONS[NUMVECTORS];
      };
    }
  }
}

#endif
