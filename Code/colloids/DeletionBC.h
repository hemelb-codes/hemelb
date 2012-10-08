// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_DELETIONBC_H
#define HEMELB_COLLOIDS_DELETIONBC_H

#include "colloids/BoundaryConditions.h"

namespace hemelb
{
  namespace colloids
  {
    class DeletionBC : public BoundaryCondition
    {
      public:
        static BoundaryCondition* ReadFromXml(io::xml::XmlAbstractionLayer& xml)
        {
          // activationDistance defaults to 1/4 lattice unit *beyond* the boundary
          LatticeDistance activationDistance = -0.25;
          xml.GetDoubleValueAndConvert("activationDistance", activationDistance);
          return new DeletionBC(activationDistance);
        }

        virtual const bool DoSomethingToParticle(
                             Particle& particle,
                             const std::vector<LatticePosition> particleToWallVectors)
        {
          // TODO: does not do *beyond* just *within* activation distance of boundary
          //LatticeDistance distance = wallNormal.GetMagnitudeSquared();
          log::Logger::Log<log::Info, log::OnePerCore>(
            "*** In DeletionBC::DoSomethingToParticle for particleId: %lu ***\n",
            particle.GetParticleId());
          return false;//distance < (activationDistance * activationDistance);
        }

        virtual const std::vector<Particle> CreateNewParticles()
        {
          return std::vector<Particle>();
        }

      protected:
        DeletionBC(LatticeDistance activationDistance) : activationDistance(activationDistance) { };

        LatticeDistance activationDistance;
    };

    class DeletionBoundaryConditionFactory : public BoundaryConditionFactory<DeletionBC> { };
  }
}
#endif /* HEMELB_COLLOIDS_DELETIONBC_H */
