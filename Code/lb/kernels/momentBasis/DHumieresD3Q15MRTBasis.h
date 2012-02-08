#ifndef HEMELB_LB_KERNELS_MOMENTBASIS_DHUMIERESD3Q15MRTBASIS_H
#define HEMELB_LB_KERNELS_MOMENTBASIS_DHUMIERESD3Q15MRTBASIS_H

#include "D3Q15.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace momentBasis
      {
        /**
         *  Class implementing the Multiple Relaxation Time (MRT) moment basis presented in in d'Humieres et al. (2002)
         *  "Multiple–relaxation–time lattice Boltzmann models in three dimensions" for the D3Q15 lattice
         */
        class DHumieresD3Q15MRTBasis
        {
          public:
            /** Moments can be separated into two groups: a) hydrodynamic (conserved) and b) kinetic (non-conserved). */
            static const unsigned NUM_KINETIC_MOMENTS = 11;

            /** Matrix used to convert from the velocities space to the reduced moment space containing only kinetic moments. */
            static const double REDUCED_MOMENT_BASIS[NUM_KINETIC_MOMENTS][D3Q15::NUMVECTORS];

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

            static void SetUpCollisionMatrix(std::vector<distribn_t>& collisionMatrix,
                                             distribn_t tau);
        };
      }
    }
  }
}
#endif //HEMELB_LB_KERNELS_MOMENTBASIS_DHUMIERESD3Q15MRTBASIS_H
