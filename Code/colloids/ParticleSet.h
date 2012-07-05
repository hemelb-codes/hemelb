#ifndef HEMELB_COLLOIDS_PARTICLESET_H
#define HEMELB_COLLOIDS_PARTICLESET_H

#include <vector>
#include "geometry/LatticeData.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
#include "mpiInclude.h"
#include "colloids/Particle.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace colloids
  {
    /** represents the set of all particles known to the local process */
    class ParticleSet
    {
      public:
        /** constructor - gets local particle information from xml config file */
        ParticleSet(const geometry::LatticeData& latDatLBM,
                    io::xml::XmlAbstractionLayer& xml,
                    lb::MacroscopicPropertyCache& propertyCache,
                    const std::vector<proc_t>& neighbourProcessors);

        /** destructor - de-allocates all Particle objects created by this Set */
        ~ParticleSet();

        /** updates the position of each particle using body forces and fluid velocity */
        const void UpdatePositions() const;

        /** calculates the effect of all body forces on each particle */
        const void CalculateBodyForces() const;

        /** calculates the effects of all particles on each lattice site */
        const void CalculateFeedbackForces() const;

        /** interpolates the fluid velocity to the location of each particle */
        const void InterpolateFluidVelocity() const;

        /** communicates the positions of all particles to&from all neighbours */
        const void CommunicateParticlePositions();

        /** communicates the partial fluid interpolations to&from all neighbours */
        const void CommunicateFluidVelocities();

      private:
        /** all particles whose position is within the local domain */
        std::vector<Particle> localParticles;

        /** temporary copy of all particles from neighbouring domains */
        std::vector<Particle> remoteParticles;

        /** contains useful geometry manipulation functions */
        const geometry::LatticeData& latDatLBM;

        /**
         * primary mechanism for interacting with the LB simulation
         * - the velocity cache: is used for velocity interpolation
         * - the forces cache  : stores the colloid feedback forces
         */
        lb::MacroscopicPropertyCache& propertyCache;

        /**
         * a vector of the processors that might be interested in
         * particles near the edge of this processor's sub-domain
         */
        const std::vector<proc_t>& neighbourProcessors;

        net::Net net;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLESET_H */
