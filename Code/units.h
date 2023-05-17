// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITS_H
#define HEMELB_UNITS_H

#include <cstdint>
#include <span>
#include "util/Vector3D.h"

namespace hemelb
{
  // Basic types for use in HemeLB.

  // ------ OLD POLICY:  ---------
  //Any variable which scales as a function of the number of sites
  // can have type site_t, processors proc_t.
  // Any variable whose precision should roughly match that of the lattice sites' velocity
  // distributions can have type distribn_t.
  typedef int64_t site_t;
  typedef int proc_t;
  typedef double distribn_t;
  typedef unsigned Direction;
  typedef uint64_t sitedata_t;

  // Span over a contiguous range of distribution-ish values
  template <std::size_t N = std::dynamic_extent>
  using ConstDistSpan = std::span<const distribn_t, N>;
  template <std::size_t N = std::dynamic_extent>
  using MutDistSpan = std::span<distribn_t, N>;

  // types used to represent domain size in blocks and for octree operations
  using U16 = std::uint16_t;
  using U64 = std::uint64_t;
  using Vec16 = util::Vector3D<U16>;

  // ------- NEW POLICY -------------
  // Types should reflect the meaning of a quantity as well as the precision
  // the type name should reflect the dimensionality and the base of the units

  typedef double PhysicalDensity;
  typedef distribn_t LatticeDensity;

  typedef double PhysicalPressure;
  typedef distribn_t LatticePressure;

  typedef double PhysicalStress;
  typedef distribn_t LatticeStress;

  typedef unsigned long LatticeTimeStep; // lattice time steps.
  typedef double LatticeTime;
  typedef double PhysicalTime; // seconds

  typedef distribn_t LatticeReciprocalTime; ///< 1/timestep
  typedef double PhysicalReciprocalTime; ///< 1/seconds
  typedef double PhysicalRate; // inverse seconds / Hz

  typedef double PhysicalMass; // kilograms

  typedef double Angle;

  typedef double PhysicalDistance; // continuous distance in physical units
  typedef double LatticeVolume; // continuous volume in physical units
  typedef double LatticeArea; // continuous area in physical units
  typedef double LatticeDistance; // continuous distance in lattice units
  typedef int64_t LatticeCoordinate; // discrete distance in lattice units

  typedef util::Vector3D<LatticeCoordinate> LatticeVector; // origin of lattice is at {0,0,0}

  // Note that Euclidean space is an affine space.
  // We can model this by distinguishing the types of points in the
  // space (positions) from the differences between them (displacements).
  // Compare to std::chrono's time_point and duration types.
  //
  // When converting between lattice and physical units, we must be aware
  // that positions are represented with respect to an origin.
  using LatticePosition = util::Vector3D<LatticeDistance>; // origin of lattice is at {0.0,0.0,0.0}
  using PhysicalPosition = util::Vector3D<PhysicalDistance>;
  using LatticeDisplacement = util::Vector3D<LatticeDistance>; // Position - Position == Displacement
  using PhysicalDisplacement = util::Vector3D<LatticeDistance>;

  typedef double PhysicalEnergy; // type for energy
  typedef double PhysicalForce; // continuous scalar force in physical units
  typedef double LatticeEnergy; // continuous scalar energy in lattice units
  typedef double LatticeForce; // continuous scalar force in lattice units
  using PhysicalModulus = double; // placeholder for any moduli, though actual dimension may differ
  typedef double LatticeModulus; // placeholder for any moduli, though actual dimension may differ
  typedef util::Vector3D<LatticeForce> LatticeForceVector; // continuous force in lattice units

  typedef double PhysicalSpeed;
  typedef double LatticeSpeed;
  typedef util::Vector3D<PhysicalSpeed> PhysicalVelocity;
  typedef util::Vector3D<LatticeSpeed> LatticeVelocity;

  using PhysicalMomentum = util::Vector3D<double>;
  using LatticeMomentum = util::Vector3D<double>;

  typedef double PhysicalPressureGradient;
  typedef double LatticePressureGradient;

  typedef double PhysicalDynamicViscosity;
  typedef double PhysicalKinematicViscosity;
  typedef double LatticeDynamicViscosity;
  typedef double LatticeKinematicViscosity;

  typedef double Dimensionless;
}
#endif //HEMELB_UNITS_H
