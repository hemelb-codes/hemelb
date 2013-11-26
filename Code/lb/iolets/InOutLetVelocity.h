//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
#include "lb/iolets/InOutLet.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetVelocity : public InOutLet
      {
        public:
          InOutLetVelocity();
          virtual ~InOutLetVelocity();
          LatticeDensity GetDensityMin() const;
          LatticeDensity GetDensityMax() const;
          LatticeDensity GetDensity(LatticeTimeStep time_step) const;

          void Reset(SimulationState &state)
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

          virtual LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const = 0;

        protected:
          LatticeDistance radius;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
