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
        UnitConverter(lb::LbmParameters* params, lb::SimulationState* state, double voxelSize);

        LatticePressure ConvertPressureToLatticeUnits(PhysicalPressure pressure) const;
        LatticeVelocity ConvertVelocityToLatticeUnits(PhysicalVelocity velocity) const;
        LatticeStress ConvertStressToLatticeUnits(PhysicalStress stress) const;
        LatticeStress ConvertPressureDifferenceToLatticeUnits(PhysicalStress pressure_grad) const;
        PhysicalPressure ConvertPressureToPhysicalUnits(LatticePressure pressure) const;
        PhysicalStress ConvertStressToPhysicalUnits(LatticeStress stress) const;
        /**
         * Templated to handle both absolute and directional velocity.
         * @param velocity
         * @return
         */
        template<class InputType>
        InputType ConvertVelocityToPhysicalUnits(InputType velocity) const
        {
          // convert velocity from lattice units to physical units (m/s)
          return velocity * latticeSpeed;
        }
        double ConvertPressureDifferenceToPhysicalUnits(distribn_t pressure_grad) const;

        PhysicalTime ConvertTimeStepToPhysicalUnits(LatticeTime time_step) const
        {
          return static_cast<double>(time_step)*simulationState->GetTimeStepLength();
        }
      private:
        lb::LbmParameters* lbmParameters;
        lb::SimulationState* simulationState;
        PhysicalLength voxelSize; //!< Lattice displacement in physical units.
        PhysicalVelocity latticeSpeed; //!< Lattice displacement length divided by time step.
    };

  }
}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
