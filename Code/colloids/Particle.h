#ifndef HEMELB_COLLOIDS_PARTICLE_H
#define HEMELB_COLLOIDS_PARTICLE_H

#include "mpiInclude.h"
#include "colloids/PersistedParticle.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace colloids
  {
    class Particle : PersistedParticle
    {
      public:
        Particle(io::xml::XmlAbstractionLayer& xml);

        util::Vector3D<double> velocity;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLE_H */
