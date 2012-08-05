// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "colloids/ParticleSet.h"
#include <algorithm>
#include "log/Logger.h"
//#include <assert.h>

namespace hemelb
{
  namespace colloids
  {
    ParticleSet::ParticleSet(const geometry::LatticeData& latDatLBM,
                             io::xml::XmlAbstractionLayer& xml,
                             lb::MacroscopicPropertyCache& propertyCache,
                             std::vector<proc_t>& neighbourProcessors) :
      localRank(topology::NetworkTopology::Instance()->GetLocalRank()),
      latDatLBM(latDatLBM),
      propertyCache(propertyCache)
    {
      std::sort(neighbourProcessors.begin(), neighbourProcessors.end());
      for (std::vector<proc_t>::const_iterator iter = neighbourProcessors.begin();
           iter != neighbourProcessors.end();
           iter++)
        scanMap.insert(scanMap.end(), scanMapContentType(*iter, scanMapElementType(0, 0)));
      scanMap.insert(scanMapContentType(localRank, scanMapElementType(0, 0)));

      // assume we are at the <Particles> node
      bool found = xml.MoveToChild("subgridParticle");
      if (found) propertyCache.velocityCache.SetRefreshFlag();
      while (found)
      {
        Particle nextParticle(latDatLBM, xml);
        if (nextParticle.isValid && nextParticle.ownerRank == localRank)
        {
          particles.push_back(nextParticle);
          scanMap[localRank].first++;
        }
        found = xml.NextSibling("subgridParticle");
      }
    };

    ParticleSet::~ParticleSet()
    {
      particles.clear();
    }

    const void ParticleSet::OutputInformation() const
    {
      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        const Particle& particle = *iter;
        particle.OutputInformation();
      }
      for (scanMapConstIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        const unsigned int& numberOfParticles = iterMap->second.first;
        const unsigned int& numberOfVelocities = iterMap->second.second;
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "ScanMap[%i] = {%i, %i}\n", neighbourRank, numberOfParticles, numberOfVelocities);
      }
    }

    const void ParticleSet::UpdatePositions()
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::ParticleSet::UpdatePositions #particles == %i ...\n",
          topology::NetworkTopology::Instance()->GetLocalRank(),
          particles.size());

      for (std::vector<Particle>::iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle& particle = *iter;
        if (particle.ownerRank == localRank)
          particle.UpdatePosition(latDatLBM);
      }
    }

    const void ParticleSet::CalculateBodyForces()
    {
      for (std::vector<Particle>::iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle& particle = *iter;
        particle.CalculateBodyForces();
      }
    }

    const void ParticleSet::CalculateFeedbackForces()
    {
      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        const Particle& particle = *iter;
        particle.CalculateFeedbackForces();
      }
    }

    const void ParticleSet::InterpolateFluidVelocity()
    {
      bool isBoring = true;
      for (std::vector<Particle>::iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle& particle = *iter;
        particle.InterpolateFluidVelocity(latDatLBM, propertyCache);
        {
          if (particle.velocity.z != 0.0)
            isBoring = false;
          if (!isBoring)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("NON-BORING velocity interpolation");
            particle.OutputInformation();
          }
        }
      }
      propertyCache.velocityCache.SetRefreshFlag();
    }

    const void ParticleSet::CommunicateParticlePositions()
    {
      /** CommunicateParticlePositions
       *    For each neighbour rank p
       *    - MPI_IRECV( number_of_remote_particles )
       *    - MPI_IRECV( list_of_remote_particles )
       *    - MPI_ISEND( number_of_local_particles )
       *    - MPI_ISEND( list_of_local_particles )
       *    MPI_WAITALL()
       *
       *  The global position of each particle is updated by the ownerRank process.
       *  The ownerRank for each particle is verified when its position is updated.
       *  Some (previously locally owned) particles may no longer be locally owned.
       *  
       */

      unsigned int& numberOfParticlesToSend = scanMap[localRank].first;
      if (scanMap.size() < 2) { return; }

      for (scanMapIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          unsigned int& numberOfParticlesToRecv = iterMap->second.first;
          net.RequestSendR(numberOfParticlesToSend, neighbourRank);
          net.RequestReceiveR(numberOfParticlesToRecv, neighbourRank);
        }
      }
      net.Dispatch();

      unsigned int numberOfParticles = 0;
      for (scanMapConstIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
        numberOfParticles += iterMap->second.first;
      particles.resize(numberOfParticles);

      std::vector<Particle>::iterator iterSendBegin = particles.begin();
      std::vector<Particle>::iterator iterRecvBegin = particles.begin() + numberOfParticlesToSend;
      for (scanMapConstIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          const unsigned int& numberOfParticlesToRecv = iterMap->second.first;
          net.RequestSend(&((PersistedParticle&)*iterSendBegin), numberOfParticlesToSend, neighbourRank);
          net.RequestReceive(&((PersistedParticle&)*(iterRecvBegin)), numberOfParticlesToRecv, neighbourRank);
          iterRecvBegin += numberOfParticlesToRecv;
        }
      }
      net.Dispatch();

      // remove particles owned by unknown ranks
      std::vector<Particle>::iterator newEndOfParticles =
        std::partition(particles.begin(), particles.end(),
                       std::bind2nd(std::mem_fun_ref(&Particle::IsOwnerRankKnown), scanMap));
      particles.erase(newEndOfParticles, particles.end());

      // sort the particles - local first, then in order of increasing owner rank
      std::sort(particles.begin(), particles.end());

      // re-build the scanMap
      for (scanMapIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
        iterMap->second.first = 0;
      for (std::vector<Particle>::const_iterator iterParticles = particles.begin();
           iterParticles != particles.end();
           iterParticles++)
        scanMap[iterParticles->ownerRank].first++;
    }

    const void ParticleSet::CommunicateFluidVelocities()
    {
      /** todo: CommunicateFluidVelocities
       *    For each neighbour rank p
       *    - MPI_IRECV( number_of_incoming_velocities )
       *    - MPI_IRECV( list_of_incoming_velocities )
       *    - MPI_ISEND( number_of_outgoing_velocities )
       *    - MPI_ISEND( list_of_outgoing_velocities )
       *    MPI_WAITALL()
       */

      if (scanMap.size() < 2) { return; }

      // exchange counts
      for (scanMapIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          unsigned int& numberOfVelocitiesToSend = iterMap->second.first;
          unsigned int& numberOfVelocitiesToRecv = iterMap->second.second;
          net.RequestSendR(numberOfVelocitiesToSend, neighbourRank);
          net.RequestReceiveR(numberOfVelocitiesToRecv, neighbourRank);
        }
      }
      net.Dispatch();

      // sum counts
      unsigned int numberOfIncomingVelocities = 0;
      for (scanMapConstIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
        numberOfIncomingVelocities += iterMap->second.second;
      velocityBuffer.resize(numberOfIncomingVelocities);

      // exchange velocities
      std::vector<Particle>::iterator iterSendBegin = particles.begin();
      std::vector<std::pair<unsigned long, util::Vector3D<double> > >::iterator
        iterRecvBegin = velocityBuffer.begin();
      for (scanMapConstIterType iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = iterMap->first;
        if (neighbourRank != localRank)
        {
          //const unsigned int& numberOfParticlesToRecv = iterMap->second.first;
          const unsigned int& numberOfVelocitiesToSend = iterMap->second.first;
          const unsigned int& numberOfVelocitiesToRecv = iterMap->second.second;
          net.RequestSend(&((Particle&)*iterSendBegin), numberOfVelocitiesToSend, neighbourRank);
          net.RequestReceive(&(*(iterRecvBegin)), numberOfVelocitiesToRecv, neighbourRank);
          iterRecvBegin += numberOfVelocitiesToRecv;
        }
      }
      net.Dispatch();

      // sum velocities
      velocityMap.clear();
      for (std::vector<std::pair<unsigned long, util::Vector3D<double> > >::const_iterator
           iterVelocityBuffer = velocityBuffer.begin();
           iterVelocityBuffer != velocityBuffer.end();
           iterVelocityBuffer++)
      {
        const unsigned long& particleId = iterVelocityBuffer->first;
        const util::Vector3D<double>& partialVelocity = iterVelocityBuffer->second;
        velocityMap[particleId] += partialVelocity;
      }

      // update particles
      for (std::vector<Particle>::iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle& particle = *iter;
        if (particle.ownerRank == localRank)
          particle.velocity += velocityMap[particle.GetParticleId()];
      }

    }

  }
}
