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
    PersistedParticle::PersistedParticle(io::xml::Element& xml, LatticeDistance voxelSize, LatticePosition geometryOrigin)
    {
      // assume we are currently at a <SubgridParticle> node
      xml.GetAttributeOrThrow("ParticleId", particleId);

      xml.GetAttributeOrThrow("InputRadiusA0", smallRadius_a0);
      // Convert to lattice units
      smallRadius_a0 = smallRadius_a0 / voxelSize;

      xml.GetAttributeOrThrow("HydrostaticRadiusAh", largeRadius_ah);
      // Convert to lattice units
      largeRadius_ah = largeRadius_ah / voxelSize;
      
      xml.GetAttributeOrThrow("Mass", mass);

      io::xml::Element initPosElem = xml.GetChildOrThrow("initialPosition");
      initPosElem.GetAttributeOrThrow("x", globalPosition.x);
      initPosElem.GetAttributeOrThrow("y", globalPosition.y);
      initPosElem.GetAttributeOrThrow("z", globalPosition.z);

      io::xml::Element fieldStrengthElem = xml.GetChildOrThrow("fieldStrength");
      fieldStrengthElem.GetAttributeOrThrow("x", fieldStrength.x);
      fieldStrengthElem.GetAttributeOrThrow("y", fieldStrength.y);
      fieldStrengthElem.GetAttributeOrThrow("z", fieldStrength.z);

      io::xml::Element temperatureElem = xml.GetChildOrThrow("temperature");
      temperatureElem.GetAttributeOrThrow("D", diffusiveTemp);
      temperatureElem.GetAttributeOrThrow("SC", softcoreTemp);

      // Convert to lattice units
      globalPosition.x = ((globalPosition.x - geometryOrigin.x) / voxelSize);
      globalPosition.y = ((globalPosition.y - geometryOrigin.y) / voxelSize);
      globalPosition.z = ((globalPosition.z - geometryOrigin.z) / voxelSize);

      lastCheckpointTimestep = 0;
      markedForDeletionTimestep = BIG_NUMBER2;
    }
    ;
  }
}
