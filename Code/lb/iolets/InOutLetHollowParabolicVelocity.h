
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETHOLLOWPARABOLICVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETHOLLOWPARABOLICVELOCITY_H
#include "lb/iolets/InOutLetVelocity.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      class InOutLetHollowParabolicVelocity : public InOutLetVelocity
      {
        public:
          InOutLetHollowParabolicVelocity();
          virtual ~InOutLetHollowParabolicVelocity();
          InOutLet* Clone() const;
          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const;

          const LatticeSpeed& GetMaxSpeed() const
          {
            return maxSpeed;
          }
          void SetMaxSpeed(const LatticeSpeed& v)
          {
            maxSpeed = v;
          }

          void SetWarmup(unsigned int warmup)
          {
            warmUpLength = warmup;
          }

          const LatticeDistance& GetInnerRadius() const
          {
            return innerRadius;
          }
          void SetInnerRadius(const LatticeDistance& r_i)
          {
            innerRadius = r_i;
          }

        protected:
          LatticeSpeed maxSpeed;
          unsigned int warmUpLength;
          LatticeDistance innerRadius;
      };
    }
  }
}
#endif // HEMELB_LB_IOLETS_INOUTLETHOLLOWPARABOLICVELOCITY_H
