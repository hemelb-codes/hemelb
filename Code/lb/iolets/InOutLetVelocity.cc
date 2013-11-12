#include "lb/iolets/InOutLetVelocity.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetVelocity::InOutLetVelocity() :
          radius(0.)
      {
      }

      InOutLetVelocity::~InOutLetVelocity()
      {
      }

      LatticeDensity InOutLetVelocity::GetDensityMin() const
      {
        return 1.0;
      }
      LatticeDensity InOutLetVelocity::GetDensityMax() const
      {
        return 1.0;
      }
      LatticeDensity InOutLetVelocity::GetDensity(LatticeTimeStep time_step) const
      {
        return 1.0;
      }
    }
  }
}
