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

    UnitConverter::UnitConverter(PhysicalTime timeStep, PhysicalDistance voxelSize,
                                 PhysicalPosition latticeOrigin) :
      latticeDistance(voxelSize), latticeTime(timeStep),
          latticeMass(BLOOD_DENSITY_Kg_per_m3 * voxelSize * voxelSize * voxelSize),
          latticeSpeed(voxelSize / latticeTime),
          latticeOrigin(latticeOrigin),
          latticePressure(latticeMass / (latticeDistance * latticeTime * latticeTime))
    {

    }

    LatticePressure UnitConverter::ConvertPressureToLatticeUnits(PhysicalPressure pressure) const
    {
      return Cs2 + ConvertPressureDifferenceToLatticeUnits(pressure - REFERENCE_PRESSURE_mmHg);
    }

    PhysicalPressure UnitConverter::ConvertPressureToPhysicalUnits(LatticePressure pressure) const
    {
      return REFERENCE_PRESSURE_mmHg + ConvertPressureDifferenceToPhysicalUnits(pressure - Cs2);
    }

    LatticePressure UnitConverter::ConvertPressureDifferenceToLatticeUnits(
                                                                           PhysicalPressure pressure_diff) const
    {
      return pressure_diff * mmHg_TO_PASCAL / latticePressure;
    }

    PhysicalPressure UnitConverter::ConvertPressureDifferenceToPhysicalUnits(
                                                                             LatticePressure pressure_diff) const
    {
      return pressure_diff * latticePressure / mmHg_TO_PASCAL;
    }

    PhysicalReciprocalTime UnitConverter::ConvertShearRateToPhysicalUnits(
                                                                          LatticeReciprocalTime shearRate) const
    {
      return shearRate / latticeTime;
    }
    LatticeDistance UnitConverter::ConvertDistanceToLatticeUnits(const PhysicalDistance& x) const
    {
      return x / latticeDistance;
    }
    PhysicalDistance UnitConverter::ConvertDistanceToPhysicalUnits(const LatticeDistance& x) const
    {
      return x * latticeDistance;
    }
    LatticePosition UnitConverter::ConvertPositionToLatticeUnits(const PhysicalPosition& x) const
    {
      return (x - latticeOrigin) / latticeDistance;
    }
    PhysicalPosition UnitConverter::ConvertPositionToPhysicalUnits(const LatticePosition& x) const
    {
      return x * latticeDistance + latticeOrigin;
    }

    LatticeSpeed UnitConverter::ConvertSpeedToLatticeUnits(const PhysicalSpeed& v) const
    {
      return v / latticeSpeed;
    }
    PhysicalSpeed UnitConverter::ConvertSpeedToPhysicalUnits(const LatticeSpeed& v) const
    {
      return v * latticeSpeed;
    }
    LatticeTime UnitConverter::ConvertTimeToLatticeUnits(const PhysicalTime& t) const
    {
      return t / latticeTime;
    }
    PhysicalTime UnitConverter::ConvertTimeToPhysicalUnits(const LatticeTime& t) const
    {
      return t * latticeTime;
    }
    PhysicalTime UnitConverter::ConvertTimeStepToPhysicalUnits(LatticeTimeStep time_step) const
    {
      return LatticeTime(time_step) * latticeTime;
    }

  }
}
