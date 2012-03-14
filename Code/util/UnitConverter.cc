#include "util/UnitConverter.h"
#include "constants.h"

namespace hemelb
{
  namespace util
  {

    UnitConverter::UnitConverter(lb::LbmParameters* iParams, lb::SimulationState* iState, double voxelSize) :
      mParams(iParams), mState(iState), voxel_size(voxelSize), CharacteristicVelocity(voxel_size/ mState->GetTimeStepLength())
    {

    }

    LatticeDensity UnitConverter::ConvertPressureToLatticeUnits(PhysicalPressure pressure) const
    {
      return Cs2 + ConvertPressureGradToLatticeUnits(pressure - REFERENCE_PRESSURE_mmHg);
    }

    PhysicalPressure UnitConverter::ConvertPressureToPhysicalUnits(LatticePressure pressure) const
    {
      return REFERENCE_PRESSURE_mmHg + ConvertPressureGradToPhysicalUnits(pressure  - Cs2);
    }

    LatticeDensity UnitConverter::ConvertPressureGradToLatticeUnits(double pressure_grad) const
    {
      return pressure_grad * mmHg_TO_PASCAL  / (CharacteristicVelocity * CharacteristicVelocity * BLOOD_DENSITY_Kg_per_m3);
    }

    double UnitConverter::ConvertPressureGradToPhysicalUnits(LatticePressure pressure_grad) const
    {
      return pressure_grad * BLOOD_DENSITY_Kg_per_m3 * CharacteristicVelocity * CharacteristicVelocity / mmHg_TO_PASCAL;
    }

    LatticeVelocity UnitConverter::ConvertVelocityToLatticeUnits(PhysicalVelocity velocity) const
    {
      return velocity / CharacteristicVelocity;
    }

    LatticeStress UnitConverter::ConvertStressToLatticeUnits(PhysicalStress stress) const
    {
      return stress / (CharacteristicVelocity*CharacteristicVelocity*BLOOD_DENSITY_Kg_per_m3);
    }

    PhysicalStress UnitConverter::ConvertStressToPhysicalUnits(PhysicalStress stress) const
    {
      // convert stress from lattice units to physical units (Pa)
      return stress * (CharacteristicVelocity*CharacteristicVelocity*BLOOD_DENSITY_Kg_per_m3);
      // stress=Force per unit area=mass * length / time^2 / length*2=mass / length * time*2
      // = mass * (length/time)^2 / length^3
    }

  }
}
