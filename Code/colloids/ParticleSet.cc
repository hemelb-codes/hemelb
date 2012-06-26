#include "colloids/ParticleSet.h"
#include "log/Logger.h"
 
namespace hemelb
{
  namespace colloids
  {
    ParticleSet::ParticleSet(const geometry::LatticeData* const latDatLBM,
                             io::xml::XmlAbstractionLayer& xml,
                             lb::MacroscopicPropertyCache& propertyCache)
    {
      // assume we are at the <Particles> node
      bool found = xml.MoveToChild("subgridParticle");
      if (found) propertyCache.velocityCache.SetRefreshFlag();
      while (found)
      {
        Particle nextParticle(latDatLBM, xml, propertyCache);
        particles.push_back(nextParticle);
        found = xml.NextSibling("subgridParticle");
      }
    };

    ParticleSet::~ParticleSet()
    {
      particles.clear();
    }

    const void ParticleSet::UpdatePositions() const
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
        log::Logger::Log<log::Debug, log::OnePerCore>(
          "[Rank: %i] In colloids::ParticleSet::UpdatePositions #particles == %i ...\n",
          topology::NetworkTopology::Instance()->GetLocalRank(),
          particles.size());

      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.UpdatePosition();
      }
    }

    const void ParticleSet::CalculateBodyForces() const
    {
      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.CalculateBodyForces();
      }
    }

    const void ParticleSet::CalculateFeedbackForces() const
    {
      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.CalculateFeedbackForces();
      }
    }

    const void ParticleSet::InterpolateFluidVelocity() const
    {
      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle particle = *iter;
        particle.InterpolateFluidVelocity();
      }
    }

    const void ParticleSet::CommunicateParticlePositions() const
    {
      /** todo: CommunicateParticlePositions
       *    For each neighbour rank p
       *    - MPI_IRECV( number_of_remote_particles )
       *    - MPI_IRECV( list_of_remote_particles )
       *    - MPI_ISEND( number_of_local_particles )
       *    - MPI_ISEND( list_of_local_particles )
       *    MPI_WAITALL()
       */
      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle particle = *iter;
      }
    }

    const void ParticleSet::CommunicateFluidVelocities() const
    {
      /** todo: CommunicateFluidVelocities
       *    For each neighbour rank p
       *    - MPI_IRECV( number_of_remote_particles )
       *    - MPI_IRECV( list_of_remote_velocities )
       *    - MPI_ISEND( number_of_local_particles )
       *    - MPI_ISEND( list_of_local_velocities )
       *    MPI_WAITALL()
       */
      for (std::vector<Particle>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle particle = *iter;
      }
    }

  }
}
