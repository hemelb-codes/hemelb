#ifndef HEMELB_COLLOIDS_PARTICLE_H
#define HEMELB_COLLOIDS_PARTICLE_H

#include "mpiInclude.h"
#include "colloids/PersistedParticle.h"
#include "geometry/LatticeData.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
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
        Particle(const geometry::LatticeData& latDatLBM,
                 io::xml::XmlAbstractionLayer& xml);

        Particle() {};

        const bool operator<(const Particle& other) const;
        const bool IsOwnerRankKnown(std::map<proc_t, unsigned int> map) const;
        const void OutputInformation() const;

        /** partial interpolation of fluid velocity - temporary value only */
        // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
        util::Vector3D<double> velocity;

        /** the effect of all body forces on this particle - this is NOT a force vector */
        // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
        util::Vector3D<double> bodyForces;

        proc_t ownerRank;
        bool isValid;

        /** updates the position of this particle using body forces and fluid velocity */
        const void UpdatePosition(const geometry::LatticeData& latDatLBM);

        /** */
        const void CalculateBodyForces();

        /** calculates the effects of all particles on each lattice site */
        const void CalculateFeedbackForces() const;

        /** interpolates the fluid velocity to the location of each particle */
        const void InterpolateFluidVelocity(
                     const geometry::LatticeData& latDatLBM,
                     const lb::MacroscopicPropertyCache& propertyCache);

        /** creates a derived MPI datatype that represents a single particle object
         *  note - this data type uses displacements rather than absolute addresses
         *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
         *  when you no longer need this type, remember to call MPI::Datatype::Free 
         */
        const MPI::Datatype CreateMpiDatatype() const;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PARTICLE_H */
