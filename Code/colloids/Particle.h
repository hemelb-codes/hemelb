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

        /** the effect of all body forces on this particle */
        util::Vector3D<double> bodyForces;

        /** updates the position of this particle using body forces and fluid velocity */
        const void UpdatePosition();

        /** */
        const void CalculateBodyForces();

        /** calculates the effects of all particles on each lattice site */
        const void CalculateFeedbackForces() const;

        /** interpolates the fluid velocity to the location of each particle */
        const void InterpolateFluidVelocity();

    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLE_H */
