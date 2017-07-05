
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_LATTICES_ADVECTIONDIFFUSIONLATTICE_H
#define HEMELB_LB_LATTICES_ADVECTIONDIFFUSIONLATTICE_H

#include "lb/lattices/Lattice.h"

namespace hemelb
{
  namespace lb
  {
    namespace lattices
    {
      template<class DmQn>
      class AdvectionDiffusionLattice : public Lattice<DmQn>
      {
        public:

          inline static void CalculateFeq(const distribn_t &concentration,
                                          const distribn_t &velocity_x,
                                          const distribn_t &velocity_y,
                                          const distribn_t &velocity_z,
                                          distribn_t f_eq[])
          {
            const distribn_t velocityMagnitudeSquared = velocity_x * velocity_x + velocity_y * velocity_y
                + velocity_z * velocity_z;

            for (Direction i = 0; i < DmQn::NUMVECTORS; ++i)
            {
              const distribn_t v_dot_ei = DmQn::CX[i] * velocity_x + DmQn::CY[i] * velocity_y
                  + DmQn::CZ[i] * velocity_z;

              f_eq[i] = DmQn::EQMWEIGHTS[i]
                  * concentration * (1 - (3. / 2.) * velocityMagnitudeSquared + (9. / 2.) * v_dot_ei * v_dot_ei
                      + 3. * v_dot_ei);
            }
          }

          inline static void CalculateConcentrationFluxFEq(const distribn_t f[],
                                                         distribn_t &concentration,
                                                         distribn_t &flux_x,
                                                         distribn_t &flux_y,
                                                         distribn_t &flux_z,
                                                         distribn_t &velocity_x,
                                                         distribn_t &velocity_y,
                                                         distribn_t &velocity_z,
                                                         distribn_t f_eq[])
          {
            Lattice<DmQn>::CalculateConcentrationAndFlux(f, concentration, flux_x, flux_y, flux_z);

            CalculateFeq(concentration, velocity_x, velocity_y, velocity_z, f_eq);
          }

      };
    }
  }
}

#endif /* HEMELB_LB_LATTICES_ADVECTIONDIFFUSIONLATTICE_H */
