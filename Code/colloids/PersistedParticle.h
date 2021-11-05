// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_PERSISTEDPARTICLE_H
#define HEMELB_COLLOIDS_PERSISTEDPARTICLE_H

#include "io/xml/XmlAbstractionLayer.h"
#include "units.h"
#include "net/MpiDataType.h"

namespace hemelb
{
  namespace colloids
  {
    /** represents a particle as stored in the xml configuration file */
    class PersistedParticle
    {
      public:
        /** constructor - gets initial values from xml configuration file */
        PersistedParticle(io::xml::Element& xml);

      protected:
        /** constructor - uses explicitly supplied values */
        PersistedParticle(unsigned long particleId, LatticeDistance a0, LatticeDistance ah,
                          PhysicalMass mass, LatticePosition globalPosition) :
            particleId(particleId), smallRadius_a0(a0), largeRadius_ah(ah), mass(mass),
                globalPosition(globalPosition)
        {
        }
        ;

        /** constructor - uses default values for each field */
        PersistedParticle()
        {
        }
        ;

        /** system-wide-unique identifier for this particle */
        unsigned long particleId;

        /** the radius of the particle */
        LatticeDistance smallRadius_a0;

        /** the hydro-static radius of the particle */
        LatticeDistance largeRadius_ah;

        /** the number of the most recent timestep during which
         *  information about this particle was written to disk
         */
        LatticeTimeStep lastCheckpointTimestep;

        /** the number of the most recent timestep during which
         *  this particle first entered the region of an outlet
         */
        LatticeTimeStep markedForDeletionTimestep;

        /** the mass of the particle */
        PhysicalMass mass;

        /** the global position of the particle in lattice units */
        LatticePosition globalPosition;
    };
  }
  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<colloids::PersistedParticle>::RegisterMpiDataType();
  }
}

#endif /* HEMELB_COLLOIDS_PERSISTEDPARTICLE_H */
