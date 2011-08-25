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
        static void Initialise(lb::LbmParameters* iParams,
                               lb::SimulationState* iState,
                               distribn_t iVoxelSize);

        static distribn_t ConvertPressureToLatticeUnits(double pressure);
        static distribn_t ConvertVelocityToLatticeUnits(double velocity);
        static distribn_t ConvertStressToLatticeUnits(double stress);
        static distribn_t ConvertPressureGradToLatticeUnits(double pressure_grad);
        static double ConvertPressureToPhysicalUnits(double pressure);
        static double ConvertStressToPhysicalUnits(double stress);
        static double ConvertVelocityToPhysicalUnits(double velocity);
        static double ConvertPressureGradToPhysicalUnits(distribn_t pressure_grad);

      private:
        static lb::LbmParameters* mParams;
        static lb::SimulationState* mState;
        static double voxel_size;
    };

  }
}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
