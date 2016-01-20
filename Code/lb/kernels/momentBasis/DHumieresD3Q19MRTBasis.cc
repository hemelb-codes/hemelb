
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cassert>
#include "lb/kernels/momentBasis/DHumieresD3Q19MRTBasis.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace momentBasis
      {

        /*
         *  Kinetic moments defined in d'Humieres et al. 2002. To get the matrix below, swap
         *  columns 9 and 11, 13 and 15, 17 and 19 to match HemeLB's lattice velocity ordering.
         *
         *  See publication for  meaning of e, epsilon, etc.
         */
        const distribn_t DHumieresD3Q19MRTBasis::REDUCED_MOMENT_BASIS[DHumieresD3Q19MRTBasis::NUM_KINETIC_MOMENTS][Lattice::NUMVECTORS] =

        { { -30, -11, -11, -11, -11, -11, -11, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 }, //e
          { 12, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, // epsilon
          { 0, -4, 4, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0 }, // q_x
          { 0, 0, 0, -4, 4, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1 }, //q_y
          { 0, 0, 0, 0, 0, -4, 4, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1 }, // q_z
          { 0, 2, 2, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2 }, // 3p_xx
          { 0, -4, -4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2 }, // 3pi_xx
          { 0, 0, 0, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0 }, // p_ww
          { 0, 0, 0, -2, -2, 2, 2, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0 }, // pi_ww
          { 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 }, // p_xy
          { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1 }, // p_yz
          { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0 }, // p_xz
          { 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0 }, // m_x
          { 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1 }, // m_y
          { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1 } //m_z
        };

        const distribn_t DHumieresD3Q19MRTBasis::BASIS_TIMES_BASIS_TRANSPOSED[NUM_KINETIC_MOMENTS] = { 2394, 252, 40,
                                                                                                       40, 40, 36, 72,
                                                                                                       12, 24, 4, 4, 4,
                                                                                                       8, 8, 8 };

        void DHumieresD3Q19MRTBasis::ProjectVelsIntoMomentSpace(const distribn_t * const velDistributions,
                                                                distribn_t * const moments)
        {
          for (unsigned momentIndex = 0; momentIndex < NUM_KINETIC_MOMENTS; momentIndex++)
          {
            moments[momentIndex] = 0.;
            for (Direction velocityIndex = 0; velocityIndex < Lattice::NUMVECTORS; velocityIndex++)
            {
              moments[momentIndex] += REDUCED_MOMENT_BASIS[momentIndex][velocityIndex]
                  * velDistributions[velocityIndex];
            }
          }
        }

        void DHumieresD3Q19MRTBasis::SetUpCollisionMatrix(std::vector<distribn_t>& collisionMatrix, distribn_t tau)
        {
          // Relaxation values taken from d'Humieres 2002, except for the kinematic viscosity where the usual tau formula is used.
          collisionMatrix.clear();
          collisionMatrix.push_back(1.19); // e (s1)
          collisionMatrix.push_back(1.4); // epsilon (s2)
          collisionMatrix.push_back(1.2); // q_x (s4)
          collisionMatrix.push_back(1.2); // q_y (s4)
          collisionMatrix.push_back(1.2); // q_z (s4)
          collisionMatrix.push_back(1.0 / tau); // 3p_xx (s9)
          collisionMatrix.push_back(1.4); // 3pi_xx s10
          collisionMatrix.push_back(1.0 / tau); // 3p_ww (s9)
          collisionMatrix.push_back(1.4); // 3pi_ww s10
          collisionMatrix.push_back(1.0 / tau); // p_xy (s13 = s9)
          collisionMatrix.push_back(1.0 / tau); // p_yz (s13 = s9)
          collisionMatrix.push_back(1.0 / tau); // p_xz (s13 = s9)
          collisionMatrix.push_back(1.98); // m_x (s16)
          collisionMatrix.push_back(1.98); // m_y (s16)
          collisionMatrix.push_back(1.98); // m_z (s16)
          assert(collisionMatrix.size() == DHumieresD3Q19MRTBasis::NUM_KINETIC_MOMENTS);
        }

      }
    }
  }
}
