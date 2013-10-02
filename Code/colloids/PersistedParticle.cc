// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
      markedForDeletionTimestep = BIG_NUMBER2;
    };
  }
}
