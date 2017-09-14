
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetHollowParabolicVelocity.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetHollowParabolicVelocity::InOutLetHollowParabolicVelocity() :
          maxSpeed(0.), warmUpLength(0), innerRadius(0.)
      {
      }

      InOutLetHollowParabolicVelocity::~InOutLetHollowParabolicVelocity()
      {
      }

      InOutLet* InOutLetHollowParabolicVelocity::Clone() const
      {
        InOutLet* copy = new InOutLetHollowParabolicVelocity(*this);
        return copy;
      }

      LatticeVelocity InOutLetHollowParabolicVelocity::GetVelocity(const LatticePosition& x,
                                                                   const LatticeTimeStep t) const
      {
        // v(r) = (-4 * vMax / ((r_o - r_i) ^ 2)) * (r - r_o) * (r - r_i)
        // where r is the distance from the centreline
        LatticePosition displ = x - position;
        LatticeDistance z = displ.Dot(normal);
        LatticeDistance r_distance = sqrt(displ.GetMagnitudeSquared() - z * z);

        // Get the max velocity
        LatticeSpeed max = maxSpeed;
        LatticeDistance r_o = radius;
        LatticeDistance r_i = innerRadius;
        // If we're in the warm-up phase, scale down the imposed velocity
        if (t < warmUpLength)
        {
          max *= t / double(warmUpLength);
        }

        // Brackets to ensure that the scalar multiplies are done before vector * scalar.
        return normal * (((-4.0 * max) / ((r_o - r_i) * (r_o - r_i))) * (r_distance - r_o) * (r_distance - r_i));
      }
    }
  }
}
