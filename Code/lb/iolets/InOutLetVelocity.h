// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
#include "lb/iolets/InOutLet.h"

namespace hemelb::lb
{
      class InOutLetVelocity : public InOutLet
      {
        public:
          InOutLetVelocity();
          ~InOutLetVelocity() override = default;
          LatticeDensity GetDensityMin() const override;
          LatticeDensity GetDensityMax() const override;
          LatticeDensity GetDensity(LatticeTimeStep time_step) const override;

          void Reset(SimulationState &state) override
          {
            //pass;
          }
          /**
           * Note that the radius and max speed for these are specified in LATTICE UNITS in the XML file.
           * This is indeed a horrible hack.
           * @return
           */
          const LatticeDistance& GetRadius() const
          {
            return radius;
          }
          void SetRadius(const LatticeDistance& r)
          {
            radius = r;
          }

          virtual LatticeVelocity GetVelocity(const LatticePosition& x,
                                              const LatticeTimeStep t) const = 0;

          //virtual LatticeVelocity GetVelocity2(const util::Vector3D<site_t> globalCoordinates,
          //                                                          const LatticeTimeStep t) const = 0;

        protected:
          LatticeDistance radius;
      };
}
#endif // HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
