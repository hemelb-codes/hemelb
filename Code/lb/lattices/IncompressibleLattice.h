
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_LATTICES_INCOMPRESSIBLELATTICE_H
#define HEMELB_LB_LATTICES_INCOMPRESSIBLELATTICE_H

#include "lb/lattices/Lattice.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<class DmQn>
      class IncompressibleLattice : public Lattice<DmQn>
      {
        public:

          inline static void CalculateFeq(const distribn_t &density,
                                          const distribn_t &momentum_x,
                                          const distribn_t &momentum_y,
                                          const distribn_t &momentum_z,
                                          distribn_t f_eq[])
          {
            const distribn_t momentumMagnitudeSquared = momentum_x * momentum_x + momentum_y * momentum_y
                + momentum_z * momentum_z;

            for (Direction i = 0; i < DmQn::NUMVECTORS; ++i)
            {
              const distribn_t mom_dot_ei = DmQn::CX[i] * momentum_x + DmQn::CY[i] * momentum_y
                  + DmQn::CZ[i] * momentum_z;

              f_eq[i] = DmQn::EQMWEIGHTS[i]
                  * (density - (3. / 2.) * momentumMagnitudeSquared + (9. / 2.) * mom_dot_ei * mom_dot_ei
                      + 3. * mom_dot_ei);
            }
          }

          inline static void CalculateDensityMomentumFEq(const distribn_t f[],
                                                         distribn_t &density,
                                                         distribn_t &momentum_x,
                                                         distribn_t &momentum_y,
                                                         distribn_t &momentum_z,
                                                         distribn_t &velocity_x,
                                                         distribn_t &velocity_y,
                                                         distribn_t &velocity_z,
                                                         distribn_t f_eq[])
          {
            Lattice<DmQn>::CalculateDensityAndMomentum(f, density, momentum_x, momentum_y, momentum_z);

            velocity_x = momentum_x;
            velocity_y = momentum_y;
            velocity_z = momentum_z;

            CalculateFeq(density, momentum_x, momentum_y, momentum_z, f_eq);
          }

          inline static bool IsLatticeCompressible()
          {
            return false;
          }

      };
    }
  }
}

#endif /* HEMELB_LB_LATTICES_INCOMPRESSIBLELATTICE_H */
