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
    /** represents the set of all particles known to the local process */
    class ParticleSet
    {
      public:
        /** constructor - gets local particle information from xml config file */
        ParticleSet(io::xml::XmlAbstractionLayer& xml);

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

      private:
        std::vector<Particle*> particles;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLESET_H */
