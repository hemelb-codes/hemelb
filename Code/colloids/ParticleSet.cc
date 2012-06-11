#include "colloids/ParticleSet.h"

namespace hemelb
{
  namespace colloids
  {
    ParticleSet::ParticleSet(io::xml::XmlAbstractionLayer& xml,
                             lb::MacroscopicPropertyCache& propertyCache)
    {
      // assume we are at the <Particles> node
      bool found = xml.MoveToChild("subgridParticle");
      while (found)
      {
        particles.push_back(new Particle(xml, propertyCache));
        found = xml.NextSibling("subgridParticle");
      }
    };

    const void ParticleSet::UpdatePositions() const
    {
      for (std::vector<Particle*>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle* particle = *iter;
        particle->UpdatePosition();
      }
    }

    const void ParticleSet::CalculateBodyForces() const
    {
      for (std::vector<Particle*>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle* particle = *iter;
        particle->CalculateBodyForces();
      }
    }

    const void ParticleSet::CalculateFeedbackForces() const
    {
      for (std::vector<Particle*>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle* particle = *iter;
        particle->CalculateFeedbackForces();
      }
    }

    const void ParticleSet::InterpolateFluidVelocity() const
    {
      for (std::vector<Particle*>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle* particle = *iter;
        particle->InterpolateFluidVelocity();
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
      for (std::vector<Particle*>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle* particle = *iter;
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
      for (std::vector<Particle*>::const_iterator iter = particles.begin();
           iter != particles.end();
           iter++)
      {
        Particle* particle = *iter;
      }
    }

  }
}
