#include "util/UnitConverter.h"
#include "constants.h"

namespace hemelb
{
  namespace util
  {

    UnitConverter::UnitConverter(lb::LbmParameters* iParams,
                                 lb::SimulationState* iState,
                                 SimConfig* iSimConfig) :
      mParams(iParams), mState(iState), voxel_size(iSimConfig->VoxelSize)
    {

    }

    distribn_t UnitConverter::ConvertPressureToLatticeUnits(double pressure) const
    {
      double temp = (PULSATILE_PERIOD_s / ((double) mState->GetTimeStepsPerCycle() * voxel_size));
      return Cs2 + (pressure - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp
          / BLOOD_DENSITY_Kg_per_m3;
    }

    double UnitConverter::ConvertPressureToPhysicalUnits(distribn_t pressure) const
    {
      double temp = ( ((double) mState->GetTimeStepsPerCycle() * voxel_size) / PULSATILE_PERIOD_s);
      return REFERENCE_PRESSURE_mmHg + ( (pressure / Cs2 - 1.0) * Cs2) * BLOOD_DENSITY_Kg_per_m3
          * temp * temp / mmHg_TO_PASCAL;
    }

    distribn_t UnitConverter::ConvertPressureGradToLatticeUnits(double pressure_grad) const
    {
      double temp = (PULSATILE_PERIOD_s / ((double) mState->GetTimeStepsPerCycle() * voxel_size));
      return pressure_grad * mmHg_TO_PASCAL * temp * temp / BLOOD_DENSITY_Kg_per_m3;
    }

    double UnitConverter::ConvertPressureGradToPhysicalUnits(distribn_t pressure_grad) const
    {
      double temp = ( ((double) mState->GetTimeStepsPerCycle() * voxel_size) / PULSATILE_PERIOD_s);
      return pressure_grad * BLOOD_DENSITY_Kg_per_m3 * temp * temp / mmHg_TO_PASCAL;
    }

    distribn_t UnitConverter::ConvertVelocityToLatticeUnits(double velocity) const
    {
      return velocity * ( ( (mParams->Tau - 0.5) / 3.0) * voxel_size) / (BLOOD_VISCOSITY_Pa_s
          / BLOOD_DENSITY_Kg_per_m3);
    }

    double UnitConverter::ConvertVelocityToPhysicalUnits(distribn_t velocity) const
    {
      // convert velocity from lattice units to physical units (m/s)
      return velocity * (BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3) / ( ( (mParams->Tau - 0.5)
          / 3.0) * voxel_size);
    }

    distribn_t UnitConverter::ConvertStressToLatticeUnits(double stress) const
    {
      return stress * (BLOOD_DENSITY_Kg_per_m3 / (BLOOD_VISCOSITY_Pa_s * BLOOD_VISCOSITY_Pa_s))
          * ( ( (mParams->Tau - 0.5) / 3.0) * voxel_size) * ( ( (mParams->Tau - 0.5) / 3.0)
          * voxel_size);
    }

    double UnitConverter::ConvertStressToPhysicalUnits(distribn_t stress) const
    {
      // convert stress from lattice units to physical units (Pa)
      return stress * BLOOD_VISCOSITY_Pa_s * BLOOD_VISCOSITY_Pa_s / (BLOOD_DENSITY_Kg_per_m3
          * ( ( (mParams->Tau - 0.5) / 3.0) * voxel_size) * ( ( (mParams->Tau - 0.5) / 3.0)
          * voxel_size));
    }

  }
}
