#ifndef HEMELB_COLLOIDS_PARTICLESET_H
#define HEMELB_COLLOIDS_PARTICLESET_H

#include <vector>
#include "mpiInclude.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "colloids/Particle.h"

namespace hemelb
{
  namespace colloids
  {
    class ParticleSet
    {
      public:
        ParticleSet(configuration::XmlAbstractionLayer& xml);
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

      private:
        std::vector<Particle*> particles;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLESET_H */
