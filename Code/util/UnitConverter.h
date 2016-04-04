
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_UNITCONVERTER_H
#define HEMELB_UTIL_UNITCONVERTER_H

#include "constants.h"
#include "units.h"
#include "util/Vector3D.h"
#include "util/Matrix3D.h"
#include "Exception.h"

namespace hemelb
{
  namespace util
  {

    class UnitConverter
    {
      public:
        UnitConverter(PhysicalTime timeStep, PhysicalDistance voxelSize, PhysicalPosition latticeOrigin);

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
          return stress * latticePressure;
        }

        /**
         * Convert a full stress tensor (including pressure and deviatoric components)
         * to physical units. Note how the diagonal is shifted by REFERENCE_PRESSURE_mmHg.
         *
         * @param stressTensor stress tensor in lattice units
         * @return stress tensor in physical units
         */
        Matrix3D ConvertFullStressTensorToPhysicalUnits(Matrix3D stressTensor) const
        {
          Matrix3D ret = stressTensor * latticePressure;
          ret.addDiagonal(REFERENCE_PRESSURE_mmHg * mmHg_TO_PASCAL);
          return ret;
        }

        /**
         * Convert a traction vector (force per unit area) to physical units. Note how a
         * REFERENCE_PRESSURE_mmHg*wallNormal component is added to account for the reference
         * pressure that was removed when converting the simulation input to lattice units.
         *
         * @param traction traction vector (computed the full stress tensor)
         * @param wallNormal wall normal at a given site
         * @return traction vector in physical units
         */
        template<class VectorType>
        Vector3D<VectorType> ConvertTractionToPhysicalUnits(Vector3D<VectorType> traction,
                                                            const Vector3D<Dimensionless>& wallNormal) const
        {
          Vector3D<VectorType> ret = traction * (latticeSpeed * latticeSpeed * BLOOD_DENSITY_Kg_per_m3);
          ret += wallNormal * REFERENCE_PRESSURE_mmHg * mmHg_TO_PASCAL;
          return ret;
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

        LatticeTime ConvertTimeToLatticeUnits(const PhysicalTime& t) const;
        PhysicalTime ConvertTimeToPhysicalUnits(const LatticeTime& t) const;
        PhysicalTime ConvertTimeStepToPhysicalUnits(LatticeTimeStep time_step) const;

        /**
         * Converts a shear rate in lattice units into physical units
         * @param shearRate shear rate in lattice units (1/time_step_length)
         * @return shear rate in physical units (1/s)
         */
        PhysicalReciprocalTime ConvertShearRateToPhysicalUnits(LatticeReciprocalTime shearRate) const;

        const PhysicalDistance& GetVoxelSize() const
        {
          return latticeDistance;
        }

        const PhysicalPosition& GetLatticeOrigin() const
        {
          return latticeOrigin;
        }
        LatticePosition GetPhysicalOrigin() const
        {
          return LatticePosition() - (latticeOrigin / latticeDistance);
        }

        template <typename T>
        T ConvertToLatticeUnits(std::string units, const T& value) const
        {
          double scale_factor;

          if (units == "m")
          {
            scale_factor = latticeDistance;
          }
          else if (units == "m/s/s")
          {
            scale_factor = latticeDistance/ (latticeTime * latticeTime);
          }
          else if (units == "N")
          {
            // F = ma so Force = mass * length / time / time and Newton = Kg * metre / second / second
            scale_factor = latticeMass * latticeDistance / (latticeTime * latticeTime);
          }
          else if (units == "s")
          {
            scale_factor = latticeTime;
          }
          else if (units == "rad" || units == "dimensionless")
          {
            scale_factor = 1.;
          }
          else if (units == "m/s")
          {
            scale_factor = latticeDistance / latticeTime;
          }
          else if (units == "mmHg")
          {
            scale_factor = latticeMass / (latticeDistance * latticeTime * latticeTime) / mmHg_TO_PASCAL;
          }
          else if (units == "Pa")
          {
            scale_factor = latticeMass / (latticeDistance * latticeTime * latticeTime);
          }
          else if (units == "mmHg/m")
          {
            scale_factor = latticeMass / (latticeDistance * latticeDistance * latticeTime * latticeTime) / mmHg_TO_PASCAL;
          }
          else if (units == "Pa/m")
          {
            scale_factor = latticeMass / (latticeDistance * latticeDistance * latticeTime * latticeTime);
          }
          else
          {
            throw Exception() << "Unknown units '" << units << "'";
          }
          return value / scale_factor;
        }

      private:
        const PhysicalDistance latticeDistance; //!< Lattice displacement in physical units.
        const PhysicalTime latticeTime;
        const PhysicalMass latticeMass;
        const PhysicalSpeed latticeSpeed; //!< Lattice displacement length divided by time step.
        const PhysicalPosition latticeOrigin;
        const PhysicalPressure latticePressure;
    };

  }
}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
