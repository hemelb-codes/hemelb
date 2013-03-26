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

    UnitConverter::UnitConverter(PhysicalTime timeStep, PhysicalDistance voxelSize, PhysicalPosition latticeOrigin) :
      voxelSize(voxelSize), timestepTime(timeStep), latticeSpeed(voxelSize / timestepTime),
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

    PhysicalReciprocalTime UnitConverter::ConvertShearRateToPhysicalUnits(LatticeReciprocalTime shearRate) const
    {
      return shearRate / timestepTime;
    }
    LatticeDistance UnitConverter::ConvertDistanceToLatticeUnits(const PhysicalDistance& x) const
    {
      return x / voxelSize;
    }
    PhysicalDistance UnitConverter::ConvertDistanceToPhysicalUnits(const LatticeDistance& x) const
    {
      return x * voxelSize;
    }
    LatticePosition UnitConverter::ConvertPositionToLatticeUnits(const PhysicalPosition& x) const
    {
      return (x - latticeOrigin) / voxelSize;
    }
    PhysicalPosition UnitConverter::ConvertPositionToPhysicalUnits(const LatticePosition& x) const
    {
      return x * voxelSize + latticeOrigin;
    }

    LatticeSpeed UnitConverter::ConvertSpeedToLatticeUnits(const PhysicalSpeed& v) const
    {
      return v / latticeSpeed;
    }
    PhysicalSpeed UnitConverter::ConvertSpeedToPhysicalUnits(const LatticeSpeed& v) const
    {
      return v * latticeSpeed;
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
