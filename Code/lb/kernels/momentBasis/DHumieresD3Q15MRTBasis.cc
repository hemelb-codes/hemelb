#include "lb/kernels/momentBasis/DHumieresD3Q15MRTBasis.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace momentBasis
      {

        /*
         *  Kinetic moments defined in d'Humieres et al. 2002. To get the matrix below, columns 8 and 9 are respectively permuted
         *  with columns 14 and 11 to match HemeLB's lattice velocity ordering.
         *
         *  See publication for  meaning of e, epsilon, etc.
         */
        const distribn_t DHumieresD3Q15MRTBasis::REDUCED_MOMENT_BASIS[DHumieresD3Q15MRTBasis::NUM_KINETIC_MOMENTS][D3Q15::NUMVECTORS] =
            { { -2, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1 }, // e
              { 16, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1 }, // epsilon
              { 0, -4, 4, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 }, // q_x
              { 0, 0, 0, -4, 4, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 }, // q_y
              { 0, 0, 0, 0, 0, -4, 4, 1, -1, -1, 1, 1, -1, -1, 1 }, // q_z
              { 0, 2, 2, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 }, // 3*p_xx
              { 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 }, // p_ww = p_yy - p_zz
              { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1 }, // p_xy
              { 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, 1, 1 }, // p_yz
              { 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1 }, // p_zx
              { 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1 } // m_xyz
            };

        const distribn_t DHumieresD3Q15MRTBasis::BASIS_TIMES_BASIS_TRANSPOSED[NUM_KINETIC_MOMENTS] =
            { 18., 360., 40., 40., 40., 12., 4., 8., 8., 8., 8. };

        void DHumieresD3Q15MRTBasis::ProjectVelsIntoMomentSpace(const distribn_t * const velDistributions,
                                                                distribn_t * const moments)
        {
          for (unsigned momentIndex = 0; momentIndex < NUM_KINETIC_MOMENTS; momentIndex++)
          {
            moments[momentIndex] = 0.;
            for (Direction velocityIndex = 0; velocityIndex < D3Q15::NUMVECTORS; velocityIndex++)
            {
              moments[momentIndex] += REDUCED_MOMENT_BASIS[momentIndex][velocityIndex]
                  * velDistributions[velocityIndex];
            }
          }
        }
      }
    }
  }
}
