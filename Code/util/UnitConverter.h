#ifndef HEMELB_UTIL_UNITCONVERTER_H
#define HEMELB_UTIL_UNITCONVERTER_H

#include "lb/LbmParameters.h"
#include "lb/SimulationState.h"
#include "constants.h"

namespace hemelb
{
  namespace util
  {

    class UnitConverter
    {
      public:
        UnitConverter(lb::LbmParameters* iParams,
                      lb::SimulationState* iState,
                      double voxelSize);

        distribn_t ConvertPressureToLatticeUnits(double pressure) const;
        distribn_t ConvertVelocityToLatticeUnits(double velocity) const;
        distribn_t ConvertStressToLatticeUnits(double stress) const;
        distribn_t ConvertPressureGradToLatticeUnits(double pressure_grad) const;
        double ConvertPressureToPhysicalUnits(double pressure) const;
        double ConvertStressToPhysicalUnits(double stress) const;
        double ConvertVelocityToPhysicalUnits(double velocity) const;
        double ConvertPressureGradToPhysicalUnits(distribn_t pressure_grad) const;

      private:
        lb::LbmParameters* mParams;
        lb::SimulationState* mState;
        double voxel_size;
    };

  }
}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
