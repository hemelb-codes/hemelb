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

#include "lb/LbmParameters.h"
#include "lb/SimulationState.h"
#include "constants.h"
#include "units.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace util
  {

    class UnitConverter
    {
      public:
        UnitConverter(lb::LbmParameters* params,
                      lb::SimulationState* state,
                      PhysicalDistance voxelSize,
                      PhysicalPosition latticeOrigin);

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
          return static_cast<double>(time_step) * simulationState->GetTimeStepLength();
        }

        /**
         * Converts a shear rate in lattice units into physical units
         * @param shearRate shear rate in lattice units (1/time_step_length)
         * @return shear rate in physical units (1/s)
         */
        PhysicalReciprocalTime ConvertShearRateToPhysicalUnits(LatticeReciprocalTime shearRate) const;

        bool Convert(std::string units, double& value) const;

        PhysicalPosition GetLatticeOrigin() const { return latticeOrigin; }
        LatticePosition GetPhysicalOrigin() const
        {
          return LatticePosition() - (latticeOrigin / voxelSize);
        }

      private:
        lb::LbmParameters* lbmParameters;
        lb::SimulationState* simulationState;
        PhysicalDistance voxelSize; //!< Lattice displacement in physical units.
        PhysicalTime timestepTime;
        PhysicalVelocity latticeSpeed; //!< Lattice displacement length divided by time step.
        PhysicalPosition latticeOrigin;
    };

  }
}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
