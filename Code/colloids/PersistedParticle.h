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
        PersistedParticle(io::xml::XmlAbstractionLayer& xml);

      protected:
        /** constructor - uses explicitly supplied values */
        PersistedParticle(unsigned long particleId, LatticeDistance a0, LatticeDistance ah,
                          LatticePosition globalPosition) :
          particleId(particleId), smallRadius_a0(a0), largeRadius_ah(ah),
          globalPosition(globalPosition)
        {};

        /** constructor - uses default values for each field */
        PersistedParticle() {};

        /** system-wide-unique identifier for this particle */
        unsigned long   particleId;

        /** the radius of the particle */
        LatticeDistance smallRadius_a0;

        /** the hydro-static radius of the particle */
        LatticeDistance largeRadius_ah;

        /** the global position of the particle in lattice units */
        LatticePosition globalPosition;

    };
  }
}

#endif /* HEMELB_COLLOIDS_PERSISTEDPARTICLE_H */
