// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UTIL_UNITCONVERTER_H
#define HEMELB_UTIL_UNITCONVERTER_H

#include "constants.h"
#include "units.h"
#include "util/Vector3D.h"
#include "util/Matrix3D.h"

namespace hemelb
{
  namespace util
  {

    class UnitConverter
    {
      public:
        UnitConverter(PhysicalTime timeStep,
                      PhysicalDistance voxelSize,
                      PhysicalPosition latticeOrigin);

        LatticePressure ConvertPressureToLatticeUnits(PhysicalPressure pressure) const;
        LatticeStress ConvertPressureDifferenceToLatticeUnits(PhysicalStress pressure_grad) const;
        PhysicalPressure ConvertPressureToPhysicalUnits(LatticePressure pressure) const;

        LatticeDistance ConvertDistanceToLatticeUnits(const PhysicalDistance& x) const;
        PhysicalDistance ConvertDistanceToPhysicalUnits(const LatticeDistance& x) const;

        LatticePosition ConvertPositionToLatticeUnits(const PhysicalPosition& x) const;
        PhysicalPosition ConvertPositionToPhysicalUnits(const LatticePosition& x) const;

        LatticeSpeed ConvertSpeedToLatticeUnits(const PhysicalSpeed& v) const;
        PhysicalSpeed ConvertSpeedToPhysicalUnits(const LatticeSpeed& v) const;

        /**
         * Convert stress from physical to lattice units, using any rank of tensor
         */
        template<class InputType>
        InputType ConvertStressToLatticeUnits(InputType stress) const
        {
          return stress / (latticeSpeed * latticeSpeed * BLOOD_DENSITY_Kg_per_m3);
        }

        /**
         * Convert stress from lattice to physical units, using any rank of tensor
         */
        template<class InputType>
        InputType ConvertStressToPhysicalUnits(InputType stress) const
        {
          return stress * (latticeSpeed * latticeSpeed * BLOOD_DENSITY_Kg_per_m3);
        }

        /**
         * Templated to handle both absolute and directional velocity.
         * @param velocity
         * @return
         */
        template<class InputType>
        InputType ConvertVelocityToLatticeUnits(InputType velocity) const
        {
          return velocity / latticeSpeed;
        }

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
          return static_cast<double>(time_step) * timestepTime;
        }

        /**
         * Converts a shear rate in lattice units into physical units
         * @param shearRate shear rate in lattice units (1/time_step_length)
         * @return shear rate in physical units (1/s)
         */
        PhysicalReciprocalTime ConvertShearRateToPhysicalUnits(LatticeReciprocalTime shearRate) const;

        bool Convert(std::string units, double& value) const;

        PhysicalPosition GetLatticeOrigin() const
        {
          return latticeOrigin;
        }
        LatticePosition GetPhysicalOrigin() const
        {
          return LatticePosition() - (latticeOrigin / voxelSize);
        }

      private:
        PhysicalDistance voxelSize; //!< Lattice displacement in physical units.
        PhysicalTime timestepTime;
        PhysicalSpeed latticeSpeed; //!< Lattice displacement length divided by time step.
        PhysicalPosition latticeOrigin;
    };

  }
}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
