#ifndef HEMELB_D3Q15_H
#define HEMELB_D3Q15_H

#include "lb/lattices/Lattice.h"

namespace hemelb
{
  class D3Q15 : public lb::lattices::Lattice<D3Q15>
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

      /** Moments can be separated into two groups: a) hydrodynamic (conserved) and b) kinetic (non-conserved). */
      static const unsigned NUM_KINETIC_MOMENTS = 11;

      /** Matrix used to convert from the velocities space to the reduced moment space containing only kinetic moments. */
      static const double REDUCED_MOMENT_BASIS[NUM_KINETIC_MOMENTS][NUMVECTORS];

      /** Diagonal matrix REDUCED_MOMENT_BASIS * REDUCED_MOMENT_BASIS'. See #61 for the MATLAB code used to compute it (in case REDUCED_MOMENT_BASIS is modified). */
      static const double BASIS_TIMES_BASIS_TRANSPOSED[NUM_KINETIC_MOMENTS];

      /**
       * Projects a velocity distributions vector into the (reduced) MRT moment space.
       *
       * @param velDistributions velocity distributions vector
       * @param moments equivalent vector in the moment space
       */
      static void ProjectVelsIntoMomentSpace(const distribn_t * const velDistributions,
                                             distribn_t * const moments);
  };
}
#endif /* HEMELB_D3Q15_H */
