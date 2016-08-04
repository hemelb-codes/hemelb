
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

	  /**
           * Calculate Feq, the orginal version
           * @param density
           * @param momentum_x
           * @param momentum_y
           * @param momentum_z
           * @param f_eq
           */
          inline static void CalculateFeq(const distribn_t &density,
                                          const distribn_t &momentum_x,
                                          const distribn_t &momentum_y,
                                          const distribn_t &momentum_z,
                                          distribn_t f_eq[])
          {
	    const distribn_t u_x = momentum_x/density;
	    const distribn_t u_y = momentum_y/density;
	    const distribn_t u_z = momentum_z/density;
	    const distribn_t u2 = u_x * u_x + u_y * u_y + u_z * u_z;

            for (Direction i = 0; i < NUMVECTORS; ++i)
            {
              const distribn_t u_dot_ei = D3Q19::CX[i] * u_x + D3Q19::CY[i] * u_y + D3Q19::CZ[i] * u_z;

              f_eq[i] = D3Q19::EQMWEIGHTS[i] * density
		* (1.
		   + 3. * u_dot_ei
		   + (9. / 2.) * u_dot_ei * u_dot_ei - (3. / 2.) * u2
		   + (27. / 2.) * D3Q19::CX[i] * u_x
		   * (  u_y * u_y * (D3Q19::CY[i] * D3Q19::CY[i] - (1./3.) )
			+ u_z * u_z * (D3Q19::CZ[i] * D3Q19::CZ[i] - (1./3.) ) )
		   + (27. / 2.) * D3Q19::CY[i] * u_y
		   * (  u_z * u_z * (D3Q19::CZ[i] * D3Q19::CZ[i] - (1./3.))
			+ u_x * u_x * (D3Q19::CX[i] * D3Q19::CX[i] - (1./3.)) )
		   + (27. / 2.) * D3Q19::CZ[i] * u_z
		   * (  u_x * u_x * (D3Q19::CX[i] * D3Q19::CX[i] - (1./3.))
			+ u_y * u_y * (D3Q19::CY[i] * D3Q19::CY[i] - (1./3.)) ) );

            }
          }

     };
    }
  }
}

#endif
