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
    /**
     * represents a single simulated biocolloid particle
     *
     * all persisted properties, i.e. those that are read in from a config file,
     * are inherited from the PersistedParticle class (which handles the I/O)
     */
    class Particle : PersistedParticle
    {
      public:
        /** constructor - gets initial values from an xml configuration file */
        Particle(io::xml::XmlAbstractionLayer& xml);

        /** partial interpolation of fluid velocity - temporary value only */
        util::Vector3D<double> velocity;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLE_H */
