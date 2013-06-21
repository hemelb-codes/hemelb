#include "lb/iolets/InOutLetVelocity.h"
#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetVelocity::InOutLetVelocity() :
        radius(0.), maxSpeed(0.), warmUpLength(0)
      {

      }

      InOutLetVelocity::~InOutLetVelocity()
      {
      }

      PhysicalPressure InOutLetVelocity::GetPressureMin() const
      {
        return REFERENCE_PRESSURE_mmHg;
      }
      PhysicalPressure InOutLetVelocity::GetPressureMax() const
      {
        return REFERENCE_PRESSURE_mmHg;
      }
      LatticeDensity InOutLetVelocity::GetDensity(LatticeTime time_step) const
      {
        return 1.0;
      }
    }
  }
}
