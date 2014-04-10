// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_PARTICLESET_H
#define HEMELB_COLLOIDS_PARTICLESET_H

#include <vector>
#include "geometry/LatticeData.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
#include "net/mpi.h"
#include "colloids/Particle.h"
#include "net/IOCommunicator.h"
#include "units.h"

namespace hemelb
{
  namespace colloids
  {
    /** represents the set of all particles known to the local process */
    class ParticleSet
    {
      public:
        /** constructor - gets local particle information from xml config file */
        ParticleSet(const geometry::LatticeData& latDatLBM,
                    io::xml::Element& xml,
                    lb::MacroscopicPropertyCache& propertyCache,
                    const hemelb::lb::LbmParameters *lbmParams,
                    std::vector<proc_t>& neighbourProcessors,
                    const net::IOCommunicator& ioComms_,
                    const std::string& outputPath);

        /** destructor - de-allocates all Particle objects created by this Set */
        ~ParticleSet();

        /** updates the position of each particle using body forces and fluid velocity */
        const void UpdatePositions();

        /** calculates the effect of all body forces on each particle */
        const void CalculateBodyForces();

        /** calculates the effects of all particles on each lattice site */
        const void CalculateFeedbackForces();

        /** applies boundary conditions to all particles **/
        const void ApplyBoundaryConditions(
                     const LatticeTimeStep currentTimestep);

        /** interpolates the fluid velocity to the location of each particle */
        const void InterpolateFluidVelocity();

        /** communicates the positions of all particles to&from all neighbours */
        const void CommunicateParticlePositions();

        /** communicates the partial fluid interpolations to&from all neighbours */
        const void CommunicateFluidVelocities();

        const void OutputInformation(const LatticeTimeStep timestep);

      private:
        const net::IOCommunicator& ioComms;
        /** cached copy of local rank (obtained from topology) */
        const proc_t localRank;

        /**
         * conatins all particles known to this process
         * they are sorted using the less than operator
         */
        std::vector<Particle> particles;

        /** map neighbourRank -> {numberOfParticlesFromThere, numberOfVelocitiesFromThere} */
        typedef std::pair<unsigned int, unsigned int> scanMapElementType;
        std::map<proc_t, scanMapElementType> scanMap;
        typedef std::map<proc_t, scanMapElementType>::const_iterator scanMapConstIterType;
        typedef std::map<proc_t, scanMapElementType>::iterator scanMapIterType;
        typedef std::pair<proc_t, scanMapElementType> scanMapContentType;
        
        /** contiguous buffer into which MPI can write all the velocities from neighbours */
        std::vector<std::pair<unsigned long, util::Vector3D<double> > > velocityBuffer;

        /** map particleId -> sumOfvelocityContributionsFromNeighbours */
        std::map<unsigned long, util::Vector3D<double> > velocityMap;

        /** contains useful geometry manipulation functions */
        const geometry::LatticeData& latDatLBM;

        /**
         * primary mechanism for interacting with the LB simulation
         * - the velocity cache  : is used for velocity interpolation
         * - the bodyForce cache : stores the colloid feedback forces
         */
        lb::MacroscopicPropertyCache& propertyCache;

        /** abstracts communication via MPI */
        net::Net net;
        /**
         * Reusable output buffer.
         */
        std::vector<char> buffer;
        /**
         * Path to write to.
         */
        const std::string& path;
        /**
         * MPI File handle to write with
         */
        net::MpiFile file;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLESET_H */
