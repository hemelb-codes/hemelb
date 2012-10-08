// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_LUBRICATIONBC_H
#define HEMELB_COLLOIDS_LUBRICATIONBC_H

#include "colloids/BoundaryConditions.h"

namespace hemelb
{
  namespace colloids
  {
    class LubricationBC : public BoundaryCondition
    {
      public:
        static BoundaryCondition* ReadFromXml(io::xml::XmlAbstractionLayer& xml)
        {
          return new LubricationBC();
        }

        virtual const bool DoSomethingToParticle(
                             Particle& particle,
                             const std::vector<LatticePosition> particleToWallVectors)
        {
          const bool keep = true;

          LatticeForceVector lubricationForce = particle.GetVelocity();
          // TODO: code lubrication formula, force = velocity.wallNormal / wallNormal.wallNormal
          // TODO: combine (add) lubrication force into particle force
          //particle.bodyForces += wallNormal;
          log::Logger::Log<log::Info, log::OnePerCore>(
            "*** In LubricationBC::DoSomethingToParticle for particleId: %lu ***\n",
            particle.GetParticleId());
          for (std::vector<LatticePosition>::const_iterator iter = particleToWallVectors.begin();
               iter != particleToWallVectors.end();
               iter++)
          {
            log::Logger::Log<log::Info, log::OnePerCore>(
              "*** In LubricationBC::DoSomethingToParticle - wall vector: {%g,%g,%g}\n",
              iter->x, iter->y, iter->z);
          }

          return keep;
        }

        virtual const std::vector<Particle> CreateNewParticles()
        {
          return std::vector<Particle>();
        }

      protected:
        LubricationBC() {};
    };

    class LubricationBoundaryConditionFactory : public BoundaryConditionFactory<LubricationBC> { };
  }
}
#endif /* HEMELB_COLLOIDS_LUBRICATIONBC_H */
