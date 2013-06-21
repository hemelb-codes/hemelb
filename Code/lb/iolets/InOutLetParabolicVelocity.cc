#include "lb/iolets/InOutLetParabolicVelocity.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetParabolicVelocity::InOutLetParabolicVelocity()
      {
      }

      InOutLetParabolicVelocity::~InOutLetParabolicVelocity()
      {
      }

      void InOutLetParabolicVelocity::DoIO(TiXmlElement *iParent, bool iIsLoading,
                                           configuration::SimConfig* iSimConfig)
      {
        iSimConfig->DoIOForParabolicVelocityInOutlet(iParent, iIsLoading, this);
      }

      InOutLet* InOutLetParabolicVelocity::Clone() const
      {
        InOutLet* copy = new InOutLetParabolicVelocity(*this);
        return copy;
      }

      LatticeVelocity InOutLetParabolicVelocity::GetVelocity(const LatticePosition& x,
                                                             const LatticeTime t) const
      {
        // v(r) = vMax (1 - r**2 / a**2)
        // where r is the distance from the centreline
        LatticePosition posLat = units->ConvertPositionToLatticeUnits(position);
        LatticePosition displ = x - posLat;
        LatticeDistance z = displ.Dot(normal);
        Dimensionless rSq = (displ.GetMagnitudeSquared() - z * z) / (radius * radius);

        // Get the max velocity
        LatticeSpeed max = maxSpeed;
        // If we're in the warm-up phase, scale down the imposed velocity
        if (t < warmUpLength)
        {
          max *= t / double(warmUpLength);
        }

        // Brackets to ensure that the scalar multiplies are done before vector * scalar.
        return normal * (max * (1. - rSq));
      }
    }
  }
}
