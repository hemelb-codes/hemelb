#ifndef HEMELB_COLLOIDS_PERSISTEDPARTICLE_H
#define HEMELB_COLLOIDS_PERSISTEDPARTICLE_H

#include "io/xml/XmlAbstractionLayer.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace colloids
  {
    /** represents a particle as stored in the xml configuration file */
    class PersistedParticle
    {
      public:
        /** constructor - gets initial values from xml configuration file */
        PersistedParticle(io::xml::XmlAbstractionLayer& xml);

      protected:
        /** constructor - uses explicitly supplied values */
        PersistedParticle(int particleId, double a0, double ah,
                          util::Vector3D<double> globalPosition) :
          particleId(particleId), smallRadius_a0(a0), largeRadius_ah(ah),
          globalPosition(globalPosition)
        {
        };

        /** system-wide-unique identifier for this particle */
        unsigned long          particleId;

        /** the radius of the particle */
        double                 smallRadius_a0;

        /** the hydro-static radius of the particle */
        double                 largeRadius_ah;

        /** the global position of the particle in lattice units */
        util::Vector3D<double> globalPosition;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PERSISTEDPARTICLE_H */
