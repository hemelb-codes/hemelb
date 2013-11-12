// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
        PersistedParticle(io::xml::Element& xml);

      protected:
        /** constructor - uses explicitly supplied values */
        PersistedParticle(unsigned long particleId,
                          LatticeDistance a0, LatticeDistance ah,
                          PhysicalMass mass,
                          LatticePosition globalPosition) :
          particleId(particleId), smallRadius_a0(a0), largeRadius_ah(ah),
          mass(mass), globalPosition(globalPosition)
        {};

        /** constructor - uses default values for each field */
        PersistedParticle() {};

        /** system-wide-unique identifier for this particle */
        unsigned long   particleId;

        /** the radius of the particle */
        LatticeDistance smallRadius_a0;

        /** the hydro-static radius of the particle */
        LatticeDistance largeRadius_ah;

        /** the number of the most recent timestep during which
         *  information about this particle was written to disk
         */
        LatticeTimeStep     lastCheckpointTimestep;

        /** the number of the most recent timestep during which
         *  this particle first entered the region of an outlet
         */
        LatticeTimeStep     markedForDeletionTimestep;

        /** the mass of the particle */
        PhysicalMass    mass;

        /** the global position of the particle in lattice units */
        LatticePosition globalPosition;
    };
  }
}

#endif /* HEMELB_COLLOIDS_PERSISTEDPARTICLE_H */
