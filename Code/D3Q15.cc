#include "D3Q15.h"
#include <cmath>

namespace hemelb
{
  const int D3Q15::CX[] = { 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 };
  const int D3Q15::CY[] = { 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 };
  const int D3Q15::CZ[] = { 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1 };
  const int* D3Q15::discreteVelocityVectors[] = { CX, CY, CZ };

  const distribn_t D3Q15::EQMWEIGHTS[] = { 2.0 / 9.0,
                                           1.0 / 9.0,
                                           1.0 / 9.0,
                                           1.0 / 9.0,
                                           1.0 / 9.0,
                                           1.0 / 9.0,
                                           1.0 / 9.0,
                                           1.0 / 72.0,
                                           1.0 / 72.0,
                                           1.0 / 72.0,
                                           1.0 / 72.0,
                                           1.0 / 72.0,
                                           1.0 / 72.0,
                                           1.0 / 72.0,
                                           1.0 / 72.0 };

  const int D3Q15::INVERSEDIRECTIONS[] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13 };
  /*
   *  Kinetic moments defined in d'Humieres 2002. To get the matrix below, columns 8 and 9 are respectively permuted
   *  with columns 14 and 11 to match HemeLB's lattice velocitiy ordering.
   *
   *  See publication for  meaning of e, epsilon, etc. Other moments definitions are possible.
   */
  const distribn_t D3Q15::REDUCED_MOMENT_BASIS[NUM_KINETIC_MOMENTS][NUMVECTORS] = {{ -2, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1}, // e
                                                                                   { 16, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1}, // epsilon
                                                                                   {  0, -4,  4,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1}, // q_x
                                                                                   {  0,  0,  0, -4,  4,  0,  0,  1, -1,  1, -1, -1,  1, -1,  1}, // q_y
                                                                                   {  0,  0,  0,  0,  0, -4,  4,  1, -1, -1,  1,  1, -1, -1,  1}, // q_z
                                                                                   {  0,  2,  2, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0}, // 3*p_xx
                                                                                   {  0,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0}, // p_ww = p_yy - p_zz
                                                                                   {  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1, -1, -1, -1, -1}, // p_xy
                                                                                   {  0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1,  1}, // p_yz
                                                                                   {  0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1}, // p_zx
                                                                                   {  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1}  // m_xyz
                                                                                  };

  const distribn_t D3Q15::BASIS_TIMES_BASIS_TRANSPOSED[NUM_KINETIC_MOMENTS] = {18., 360., 40., 40., 40., 12., 4., 8., 8., 8., 8.};



  void D3Q15::ProjectVelsIntoMomentSpace(const distribn_t * const velDistributions,
                                         distribn_t * const moments)
  {
    for (unsigned momentIndex = 0; momentIndex < NUM_KINETIC_MOMENTS; momentIndex++)
    {
      moments[momentIndex] = 0.;
      for (unsigned velocityIndex = 0; velocityIndex < NUMVECTORS; velocityIndex++)
      {
        moments[momentIndex] += REDUCED_MOMENT_BASIS[momentIndex][velocityIndex]
            * velDistributions[velocityIndex];
      }
    }
  }

}
