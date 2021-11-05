// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/PersistedParticle.h"
#include "constants.h"

namespace hemelb
{
  namespace colloids
  {
    PersistedParticle::PersistedParticle(io::xml::Element& xml)
    {
      // assume we are currently at a <SubgridParticle> node
      xml.GetAttributeOrThrow("ParticleId", particleId);
      // TODO: Convert to lattice
      xml.GetAttributeOrThrow("InputRadiusA0", smallRadius_a0);
      // TODO: Convert to lattice
      xml.GetAttributeOrThrow("HydrostaticRadiusAh", largeRadius_ah);
      // TODO: This MAY not need to be converted to lattice (wasn't before my refactor)
      xml.GetAttributeOrThrow("Mass", mass);
      // TODO: Convert to lattice
      io::xml::Element initPosElem = xml.GetChildOrThrow("initialPosition");
      initPosElem.GetAttributeOrThrow("x", globalPosition.x);
      initPosElem.GetAttributeOrThrow("y", globalPosition.y);
      initPosElem.GetAttributeOrThrow("z", globalPosition.z);

      lastCheckpointTimestep = 0;
      markedForDeletionTimestep = SITE_OR_BLOCK_SOLID;
    };
  }
}
