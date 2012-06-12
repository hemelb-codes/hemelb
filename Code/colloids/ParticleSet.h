#ifndef HEMELB_COLLOIDS_PARTICLESET_H
#define HEMELB_COLLOIDS_PARTICLESET_H

#include <vector>
#include "mpiInclude.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
#include "colloids/Particle.h"

namespace hemelb
{
  namespace colloids
  {
    /** represents the set of all particles known to the local process */
    class ParticleSet
    {
      public:
        /** constructor - gets local particle information from xml config file */
        ParticleSet(io::xml::XmlAbstractionLayer& xml,
                    lb::MacroscopicPropertyCache& propertyCache);

        /** destructor - de-allocates all Particle objects created by this Set */
        ~ParticleSet()
        {
          for (std::vector<Particle*>::const_iterator iter = particles.begin();
               iter != particles.end();
               iter++)
          {
            delete *iter;
          }
          particles.clear();
        }

        /** updates the position of each particle using body forces and fluid velocity */
        const void UpdatePositions() const;

        /** calculates the effect of all body forces on each particle */
        const void CalculateBodyForces() const;

        /** calculates the effects of all particles on each lattice site */
        const void CalculateFeedbackForces() const;

        /** interpolates the fluid velocity to the location of each particle */
        const void InterpolateFluidVelocity() const;

        /** communicates the positions of all particles to&from all neighbours */
        const void CommunicateParticlePositions() const;

        /** communicates the partial fluid interpolations to&from all neighbours */
        const void CommunicateFluidVelocities() const;

      private:
        std::vector<Particle*> particles;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLESET_H */
