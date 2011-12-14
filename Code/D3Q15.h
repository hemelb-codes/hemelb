#ifndef HEMELB_D3Q15_H
#define HEMELB_D3Q15_H

#include "lb/lattices/Lattice.h"

namespace hemelb
{
  class D3Q15 : public lb::lattices::Lattice<D3Q15>
  {
    public:
      // The number of discrete velocity vectors
      static const unsigned int NUMVECTORS = 15;

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
#endif /* HEMELB_D3Q15_H */
