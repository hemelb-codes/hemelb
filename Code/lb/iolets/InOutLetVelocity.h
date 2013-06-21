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
          PhysicalPressure GetPressureMin() const;
          PhysicalPressure GetPressureMax() const;
          LatticeDensity GetDensity(LatticeTime time_step) const;

          void Reset(SimulationState &state)
          {
            //pass;
          }
          /**
           * Note that the radius and max speed for these are specified in LATTICE UNITS in the XML file.
           * This is indeed a horrible hack.
           * @return
           */
          LatticeDistance& GetRadius()
          {
            return radius;
          }
          void SetRadius(LatticeDistance r)
          {
            radius = r;
          }

          LatticeSpeed& GetMaxSpeed()
          {
            return maxSpeed;
          }

          void SetMaxSpeed(LatticeSpeed v)
          {
            maxSpeed = v;
          }

          void SetWarmup(unsigned int warmup)
          {
            warmUpLength = warmup;
          }

          virtual LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTime t) const = 0;

        protected:
          LatticeDistance radius;
          LatticeSpeed maxSpeed;
          unsigned int warmUpLength;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
