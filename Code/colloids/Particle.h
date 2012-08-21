// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_PARTICLE_H
#define HEMELB_COLLOIDS_PARTICLE_H

#include "mpiInclude.h"
#include "colloids/PersistedParticle.h"
#include "geometry/LatticeData.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace colloids
  {
    /**
     * represents a single simulated biocolloid particle
     *
     * all persisted properties, i.e. those that are read in from a config file,
     * are inherited from the PersistedParticle class (which handles the I/O)
     */
    class Particle : PersistedParticle
    {
      public:
        /** constructor - gets initial values from an xml configuration file */
        Particle(const geometry::LatticeData& latDatLBM,
                 io::xml::XmlAbstractionLayer& xml);

        /** constructor - gets an invalid particle for making MPI data types */
        Particle() {};

        /** property getter for particleId */
        const unsigned long GetParticleId() const { return particleId; }
        const LatticePosition& GetGlobalPosition() const { return globalPosition; }
        const PhysicalMass GetMass() const {return mass; }

        /** property getter for ownerRank */
        const proc_t GetOwnerRank() const { return ownerRank; }

        /** property getter for isValid */
        const bool IsValid() const { return isValid; }

        /**
         * less than operator for comparing particle objects
         *
         * when used to sort a container of particle objects
         * the ordering produced by this operator is:
         * - increasing particleId
         * - grouped by owner rank
         * - with local rank first
         */
        const bool operator<(const Particle& other) const;

        /** determines if the owner rank of this particle is an existing key in map */
        const bool IsOwnerRankKnown(std::map<proc_t, std::pair<unsigned int, unsigned int> > map) const;

        /** for debug purposes only - outputs all properties to info log */
        const void OutputInformation() const;

        /** obtains the fluid viscosity at the position of this particle */
        // TODO: currently returns BLOOD_VISCOSITY_Pa_s, which has the wrong units
        const DimensionlessQuantity GetViscosity() const;

        /** calculates the drag coefficient = 1/(6*pi*viscosity*radius) */
        const DimensionlessQuantity CalculateDrag();

        /** updates the position of this particle using body forces and fluid velocity */
        const void UpdatePosition(const geometry::LatticeData& latDatLBM);

        /** calculates the effects of all body forces on this particle */
        const void CalculateBodyForces();

        /** calculates the effects of this particle on each lattice site */
        const void CalculateFeedbackForces(const geometry::LatticeData& latDatLBM) const;

        /** interpolates the fluid velocity to the location of each particle */
        const void InterpolateFluidVelocity(
                     const geometry::LatticeData& latDatLBM,
                     const lb::MacroscopicPropertyCache& propertyCache);

        /** accumulate contributions to velocity from remote processes */
        const void AccumulateVelocity(util::Vector3D<double>& contribution)
        {
          velocity += contribution;
        };

        /** creates a derived MPI datatype that represents a single particle object
         *  the fields included are all those from the PersistedParticle base class
         *  note - this data type uses displacements rather than absolute addresses
         *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
         *  when you no longer need this type, remember to call MPI::Datatype::Free 
         */
        const MPI::Datatype CreateMpiDatatypeWithPosition() const;

        /** creates a derived MPI datatype that represents a single particle object
         *  the fields included in this type are: particleId and velocity(xyz) only
         *  note - this data type uses displacements rather than absolute addresses
         *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
         *  when you no longer need this type, remember to call MPI::Datatype::Free 
         */
        const MPI::Datatype CreateMpiDatatypeWithVelocity() const;

      private:
        /** partial interpolation of fluid velocity - temporary value only */
        // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
        util::Vector3D<double> velocity;

        /** the effect of all body forces on this particle - this is NOT a force vector */
        // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
        util::Vector3D<double> bodyForces;

        DimensionlessQuantity drag;

        proc_t ownerRank;

        bool isValid;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLE_H */
