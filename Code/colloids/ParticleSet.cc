#include "colloids/ParticleSet.h"
#include <algorithm>
#include "log/Logger.h"
#undef NDEBUG
#include <assert.h>

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
      //, neighbourProcessors(neighbourProcessors)
    {
      std::sort(neighbourProcessors.begin(), neighbourProcessors.end());
      for (std::vector<proc_t>::const_iterator iter = neighbourProcessors.begin();
           iter != neighbourProcessors.end();
           iter++)
        scanMap.insert(scanMap.end(), std::pair<proc_t, unsigned int>(*iter, 0));
      scanMap.insert(std::pair<proc_t, unsigned int>(localRank, 0));

      // assume we are at the <Particles> node
      bool found = xml.MoveToChild("subgridParticle");
      if (found) propertyCache.velocityCache.SetRefreshFlag();
      while (found)
      {
        Particle nextParticle(latDatLBM, xml);
        //localParticles.push_back(nextParticle);
        if (nextParticle.isValid && nextParticle.ownerRank == localRank)
        {
          particles.push_back(nextParticle);
          scanMap[localRank]++;
        }
        found = xml.NextSibling("subgridParticle");
      }
    };

    ParticleSet::~ParticleSet()
    {
      //localParticles.clear();
      //remoteParticles.clear();
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
      std::map<proc_t, unsigned int>::const_iterator iterMap = scanMap.begin();
      for (iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = (*iterMap).first;
        const unsigned int& numberOfParticles = (*iterMap).second;
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "ScanMap[%i] = %i\n", neighbourRank, numberOfParticles);
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
        if (particle.velocity.z != 0.0)
          isBoring = false;
        if (!isBoring)
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("NON-BORING velocity interpolation");
          particle.OutputInformation();
        }
      }
      propertyCache.velocityCache.SetRefreshFlag();
    }

    const void ParticleSet::CommunicateParticlePositions()
    {
      /** todo: CommunicateParticlePositions
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

      std::map<proc_t, unsigned int>::iterator iterMap = scanMap.begin();
      unsigned int& numberOfParticlesToSend = scanMap[localRank];
      if (scanMap.size() < 2) { return; }

      for (iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = (*iterMap).first;
        if (neighbourRank != localRank)
        {
          unsigned int& numberOfParticlesToRecv = (*iterMap).second;
          net.RequestSendR(numberOfParticlesToSend, neighbourRank);
          net.RequestReceiveR(numberOfParticlesToRecv, neighbourRank);
        }
      }
      net.Dispatch();

      unsigned int numberOfParticles = 0;
      for (iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
        numberOfParticles += (*iterMap).second;
      particles.resize(numberOfParticles);

      std::vector<Particle>::iterator iterSendBegin = particles.begin();
      std::vector<Particle>::iterator iterRecvBegin = particles.begin() + numberOfParticlesToSend;
      for (iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
      {
        const proc_t& neighbourRank = (*iterMap).first;
        if (neighbourRank != localRank)
        {
          const unsigned int& numberOfParticlesToRecv = (*iterMap).second;
          net.RequestSend(&(*iterSendBegin), numberOfParticlesToSend, neighbourRank);
          net.RequestReceive(&(*(iterRecvBegin)), numberOfParticlesToRecv, neighbourRank);
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
      for (iterMap = scanMap.begin();
           iterMap != scanMap.end();
           iterMap++)
        (*iterMap).second = 0;
      for (std::vector<Particle>::const_iterator iterParticles = particles.begin();
           iterParticles != particles.end();
           iterParticles++)
      {
        scanMap[(*iterParticles).ownerRank]++;
      }

/*
      unsigned int numberOfParticlesToSend = localParticles.size();
      unsigned int numberOfParticlesToReceive[55];//neighbourProcessors.size() + 1];
      numberOfParticlesToReceive[0] = 0;

      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::ParticleSet::CommunicateParticlePositions #particles == %i, #neighbours == %i ...\n",
          numberOfParticlesToSend, 55);//neighbourProcessors.size());

      //if (neighbourProcessors.empty()) { return; }

      unsigned int iterNeighbourIndex = 0;
      for (std::vector<proc_t>::const_iterator iterNeighbourRank = neighbourProcessors.begin();
           iterNeighbourRank != neighbourProcessors.end();
           iterNeighbourRank++)
      {
        if (log::Logger::ShouldDisplay<log::Debug>())
          log::Logger::Log<log::Debug, log::OnePerCore>(
            "In colloids::ParticleSet::CommunicateParticlePositions sending (int) %i to rank %i\n",
            numberOfParticlesToSend, *iterNeighbourRank);

        net.RequestSendR(numberOfParticlesToSend, *iterNeighbourRank);
        iterNeighbourIndex++;
        net.RequestReceive(&numberOfParticlesToReceive[iterNeighbourIndex], 1, *iterNeighbourRank);

        if (log::Logger::ShouldDisplay<log::Debug>())
          log::Logger::Log<log::Debug, log::OnePerCore>(
            "In colloids::ParticleSet::CommunicateParticlePositions recving (int) [%i,1] from rank %i\n",
            iterNeighbourIndex, *iterNeighbourRank);
      }
      net.Dispatch();

      unsigned int numberOfRemoteParticles = 0;
      for (iterNeighbourIndex = 0;
           iterNeighbourIndex < neighbourProcessors.size();
           iterNeighbourIndex++)
      {
        numberOfRemoteParticles += numberOfParticlesToReceive[iterNeighbourIndex + 1];
      }
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::ParticleSet::CommunicateParticlePositions #remoteParticles == %i ...\n",
          numberOfRemoteParticles);

      if (remoteParticles.size() < numberOfRemoteParticles)
        remoteParticles.resize(numberOfRemoteParticles);

      iterNeighbourIndex = 0;
      for (std::vector<proc_t>::const_iterator iterNeighbourRank = neighbourProcessors.begin();
           iterNeighbourRank != neighbourProcessors.end();
           iterNeighbourRank++)
      {
        if (log::Logger::ShouldDisplay<log::Debug>())
          log::Logger::Log<log::Debug, log::OnePerCore>(
            "In colloids::ParticleSet::CommunicateParticlePositions sending (Particle) [%i] to rank %i\n",
            localParticles.size(), *iterNeighbourRank);

        net.RequestSendV(localParticles, *iterNeighbourRank);
        net.RequestReceive(&(remoteParticles[numberOfParticlesToReceive[iterNeighbourIndex]]),
                           numberOfParticlesToReceive[iterNeighbourIndex + 1],
                           *iterNeighbourRank);

        if (log::Logger::ShouldDisplay<log::Debug>())
          log::Logger::Log<log::Debug, log::OnePerCore>(
            "In colloids::ParticleSet::CommunicateParticlePositions recving (Particle) [%i,%i] from rank %i\n",
            numberOfParticlesToReceive[iterNeighbourIndex],
            numberOfParticlesToReceive[iterNeighbourIndex + 1],
            *iterNeighbourRank);

        numberOfParticlesToReceive[iterNeighbourIndex + 1] +=
          numberOfParticlesToReceive[iterNeighbourIndex];
        iterNeighbourIndex++;
      }
      net.Dispatch();

      std::sort(remoteParticles.begin(), remoteParticles.end());
*/
    }

    const void ParticleSet::CommunicateFluidVelocities()
    {
      /** todo: CommunicateFluidVelocities
       *    For each neighbour rank p
       *    - MPI_IRECV( number_of_remote_particles )
       *    - MPI_IRECV( list_of_remote_velocities )
       *    - MPI_ISEND( number_of_local_particles )
       *    - MPI_ISEND( list_of_local_velocities )
       *    MPI_WAITALL()
       */
      for (std::vector<Particle>::iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle& particle = *iter;
      }
    }

  }
}
