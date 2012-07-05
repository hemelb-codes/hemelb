#ifndef HEMELB_COLLOIDS_PERSISTEDPARTICLE_H
#define HEMELB_COLLOIDS_PERSISTEDPARTICLE_H

#include "io/xml/XmlAbstractionLayer.h"
#include "units.h"

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
        PersistedParticle(unsigned long particleId, LatticeDistance a0, LatticeDistance ah,
                          LatticePosition globalPosition) :
          particleId(particleId), smallRadius_a0(a0), largeRadius_ah(ah),
          globalPosition(globalPosition)
        {
        };

        /** system-wide-unique identifier for this particle */
        unsigned long   particleId;

        /** the radius of the particle */
        LatticeDistance smallRadius_a0;

        /** the hydro-static radius of the particle */
        LatticeDistance largeRadius_ah;

        /** the global position of the particle in lattice units */
        LatticePosition globalPosition;

        PersistedParticle() {};
    };
  }
}

#endif /* HEMELB_COLLOIDS_PERSISTEDPARTICLE_H */
