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
        UnitConverter(lb::LbmParameters* iParams, lb::SimulationState* iState, double voxelSize);

        distribn_t ConvertPressureToLatticeUnits(double pressure) const;
        double ConvertVelocityToLatticeUnits(double velocity) const;
        distribn_t ConvertStressToLatticeUnits(double stress) const;
        distribn_t ConvertPressureGradToLatticeUnits(double pressure_grad) const;
        double ConvertPressureToPhysicalUnits(double pressure) const;
        double ConvertStressToPhysicalUnits(double stress) const;
        /**
         * Templated to handle both absolute and directional velocity.
         * @param velocity
         * @return
         */
        template<class InputType>
        InputType ConvertVelocityToPhysicalUnits(InputType velocity) const
        {
          // convert velocity from lattice units to physical units (m/s)
          return velocity * (BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3) / ( ( (mParams->GetTau() - 0.5) / 3.0)
              * voxel_size);
        }
        double ConvertPressureGradToPhysicalUnits(distribn_t pressure_grad) const;

        PhysicalTime ConvertTimeStepToPhysicalUnits(LatticeTime time_step) const
        {
          return static_cast<double>(time_step)*mState->GetTimeStepLength();
        }
      private:
        lb::LbmParameters* mParams;
        lb::SimulationState* mState;
        double voxel_size;
    };

  }
}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
