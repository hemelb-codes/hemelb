// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "util/UnitConverter.h"
#include "constants.h"

namespace hemelb
{
  namespace util
  {

    UnitConverter::UnitConverter(lb::LbmParameters* params,
                                 lb::SimulationState* state,
                                 PhysicalDistance voxelSize,
                                 PhysicalPosition latticeOrigin) :
      lbmParameters(params), simulationState(state),
      voxelSize(voxelSize), timestepTime(simulationState->GetTimeStepLength()),
      latticeSpeed(voxelSize / timestepTime),
      latticeOrigin(latticeOrigin)
    {

    }

    LatticeDensity UnitConverter::ConvertPressureToLatticeUnits(PhysicalPressure pressure) const
    {
      return Cs2 + ConvertPressureDifferenceToLatticeUnits(pressure - REFERENCE_PRESSURE_mmHg);
    }

    PhysicalPressure UnitConverter::ConvertPressureToPhysicalUnits(LatticePressure pressure) const
    {
      return REFERENCE_PRESSURE_mmHg + ConvertPressureDifferenceToPhysicalUnits(pressure - Cs2);
    }

    LatticeDensity UnitConverter::ConvertPressureDifferenceToLatticeUnits(double pressure_grad) const
    {
      return pressure_grad * mmHg_TO_PASCAL / (latticeSpeed * latticeSpeed * BLOOD_DENSITY_Kg_per_m3);
    }

    double UnitConverter::ConvertPressureDifferenceToPhysicalUnits(LatticePressure pressure_grad) const
    {
      return pressure_grad * BLOOD_DENSITY_Kg_per_m3 * latticeSpeed * latticeSpeed / mmHg_TO_PASCAL;
    }

    LatticeSpeed UnitConverter::ConvertSpeedToLatticeUnits(PhysicalSpeed speed) const
    {
      return speed / latticeSpeed;
    }

    LatticeStress UnitConverter::ConvertStressToLatticeUnits(PhysicalStress stress) const
    {
      return stress / (latticeSpeed * latticeSpeed * BLOOD_DENSITY_Kg_per_m3);
    }

    PhysicalStress UnitConverter::ConvertStressToPhysicalUnits(PhysicalStress stress) const
    {
      // convert stress from lattice units to physical units (Pa)
      return stress * (latticeSpeed * latticeSpeed * BLOOD_DENSITY_Kg_per_m3);
      // stress=Force per unit area=mass * length / time^2 / length^2=mass / length * time^2
      // = mass * (length/time)^2 / length^3
    }

    PhysicalReciprocalTime UnitConverter::ConvertShearRateToPhysicalUnits(LatticeReciprocalTime shearRate) const
    {
      return shearRate / simulationState->GetTimeStepLength();
    }

    bool UnitConverter::Convert(std::string units, double& value) const
    {
      if (units == "metres" || units == "m")
        value /= voxelSize;
      else if (units == "si_acceleration" || units == "m/s/s" || units == "ms-2")
        value *= timestepTime * timestepTime / voxelSize;
      else if (units == "Newtons" || units == "N")
        // F = ma so Force = mass * length / time / time and Newton = Kg * metre / second / second
        value *= voxelSize * voxelSize * latticeSpeed * latticeSpeed * BLOOD_DENSITY_Kg_per_m3;
      else
        return false;
      return true;
    }
  }
}
