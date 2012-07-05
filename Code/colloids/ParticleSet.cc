#include "colloids/ParticleSet.h"
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
                             const std::vector<proc_t>& neighbourProcessors) :
      latDatLBM(latDatLBM),
      propertyCache(propertyCache),
      neighbourProcessors(neighbourProcessors)
    {
      // assume we are at the <Particles> node
      bool found = xml.MoveToChild("subgridParticle");
      if (found) propertyCache.velocityCache.SetRefreshFlag();
      while (found)
      {
        Particle nextParticle(latDatLBM, xml);
        localParticles.push_back(nextParticle);
        found = xml.NextSibling("subgridParticle");
      }
    };

    ParticleSet::~ParticleSet()
    {
      localParticles.clear();
      remoteParticles.clear();
    }

    const void ParticleSet::UpdatePositions() const
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "[Rank: %i] In colloids::ParticleSet::UpdatePositions #particles == %i ...\n",
          topology::NetworkTopology::Instance()->GetLocalRank(),
          localParticles.size());

      for (std::vector<Particle>::const_iterator iter = localParticles.begin();
           iter != localParticles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.UpdatePosition();
      }
    }

    const void ParticleSet::CalculateBodyForces() const
    {
      for (std::vector<Particle>::const_iterator iter = localParticles.begin();
           iter != localParticles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.CalculateBodyForces();
      }
    }

    const void ParticleSet::CalculateFeedbackForces() const
    {
      for (std::vector<Particle>::const_iterator iter = localParticles.begin();
           iter != localParticles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.CalculateFeedbackForces();
      }
    }

    const void ParticleSet::InterpolateFluidVelocity() const
    {
      for (std::vector<Particle>::const_iterator iter = localParticles.begin();
           iter != localParticles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.InterpolateFluidVelocity(latDatLBM, propertyCache);
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
       */

      unsigned int numberOfParticlesToSend = localParticles.size();
      unsigned int numberOfParticlesToReceive[neighbourProcessors.size() + 1];
      numberOfParticlesToReceive[0] = 0;

      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "In colloids::ParticleSet::CommunicateParticlePositions #particles == %i, #neighbours == %i ...\n",
          numberOfParticlesToSend, neighbourProcessors.size());

      //assert(neighbourProcessors.size());
      if (neighbourProcessors.empty()) { return; }

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

      {
        Particle temp[2];
        unsigned int extentExpected = &(temp[1])-&(temp[0]);
        for (int iterNeighbourIndex=0;
           iterNeighbourIndex < neighbourProcessors.size();
           iterNeighbourIndex++)
        {
          unsigned int extentActual = &(remoteParticles[iterNeighbourIndex + 1])
                                    - &(remoteParticles[iterNeighbourIndex]);
          assert(extentActual == extentExpected);
        }
      }

      iterNeighbourIndex = 0;
      for (std::vector<proc_t>::const_iterator iterNeighbourRank = neighbourProcessors.begin();
           iterNeighbourRank != neighbourProcessors.end();
           iterNeighbourRank++)
      {
        if (log::Logger::ShouldDisplay<log::Debug>())
          log::Logger::Log<log::Debug, log::OnePerCore>(
            "In colloids::ParticleSet::CommunicateParticlePositions sending (Particle) %i to rank %i\n",
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

      assert(remoteParticles.size() == numberOfRemoteParticles);
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
      for (std::vector<Particle>::const_iterator iter = localParticles.begin();
           iter != localParticles.end();
           iter++)
      {
        Particle particle = *iter;
      }
    }

  }
}
