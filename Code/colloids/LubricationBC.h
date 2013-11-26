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
#include "io/xml/XmlAbstractionLayer.h"

namespace hemelb
{
  namespace colloids
  {
    class LubricationBC : public BoundaryCondition
    {
      public:
        static BoundaryCondition* ReadFromXml(io::xml::Element& xml)
        {
          LatticeDistance effectiveRange = 1.0;
          xml.GetAttributeOrThrow("effectiveRange", effectiveRange);
          return new LubricationBC(effectiveRange);
        }

        virtual const bool DoSomethingToParticle(
                             Particle& particle,
                             const std::vector<LatticePosition> particleToWallVectors)
        {
          const bool keep = true;

          log::Logger::Log<log::Trace, log::OnePerCore>(
            "*** In LubricationBC::DoSomethingToParticle for particleId: %lu ***\n",
            particle.GetParticleId());

          LatticeVelocity lubricationVelocityAdjustment;
          particle.SetLubricationVelocityAdjustment(lubricationVelocityAdjustment);

          for (std::vector<LatticePosition>::const_iterator iter = particleToWallVectors.begin();
               iter != particleToWallVectors.end();
               iter++)
          {
            const LatticePosition particleToWallVector = *iter;
            const LatticeDistance separation_h = particleToWallVector.GetMagnitude()
                                               - particle.GetRadius();

            log::Logger::Log<log::Trace, log::OnePerCore>(
              "*** In LubricationBC::DoSomethingToParticle - wall vector: {%g,%g,%g}, mag: %g, particle radius: %g, separation_h: %g\n",
              particleToWallVector.x,
              particleToWallVector.y,
              particleToWallVector.z,
              particleToWallVector.GetMagnitude(),
              particle.GetRadius(),
              separation_h);

            if (separation_h <= effectiveRange)
            {
              lubricationVelocityAdjustment +=
                particleToWallVector.GetNormalised()
                * particle.GetVelocity().Dot(particleToWallVector.GetNormalised())
                * ( (separation_h - effectiveRange) / (separation_h * effectiveRange) )
                * particle.GetRadius() * particle.GetRadius()
                * particle.GetInverseNormalisedRadius();

              log::Logger::Log<log::Trace, log::OnePerCore>(
                "*** In LubricationBC::DoSomethingToParticle - radius: %g, separation: %g, adj: {%g,%g,%g}\n",
                particle.GetInverseNormalisedRadius(),
                ( (effectiveRange - separation_h) / (separation_h * effectiveRange) ),
                lubricationVelocityAdjustment.x,
                lubricationVelocityAdjustment.y,
                lubricationVelocityAdjustment.z);
            } else {
              log::Logger::Log<log::Trace, log::OnePerCore>(
                "*** In LubricationBC::DoSomethingToParticle - separation: %g, range: %g\n",
                separation_h, effectiveRange);
            }
          }
          particle.SetLubricationVelocityAdjustment(lubricationVelocityAdjustment);

          log::Logger::Log<log::Trace, log::OnePerCore>(
            "*** In LubricationBC::DoSomethingToParticle - particleId: %lu, vel before: {%g,%g,%g}, total adj: {%g,%g,%g}, vel after: {%g,%g,%g}\n",
            particle.GetParticleId(),
            particle.GetVelocity().x - lubricationVelocityAdjustment.x,
            particle.GetVelocity().y - lubricationVelocityAdjustment.y,
            particle.GetVelocity().z - lubricationVelocityAdjustment.z,
            lubricationVelocityAdjustment.x,
            lubricationVelocityAdjustment.y,
            lubricationVelocityAdjustment.z,
            particle.GetVelocity().x,
            particle.GetVelocity().y,
            particle.GetVelocity().z);

          return keep;
        }

        virtual const std::vector<Particle> CreateNewParticles()
        {
          return std::vector<Particle>();
        }

      protected:
        LubricationBC(const LatticeDistance effectiveRange) : effectiveRange(effectiveRange) {};
        LatticeDistance effectiveRange;
    };

    class LubricationBoundaryConditionFactory : public BoundaryConditionFactory<LubricationBC> { };
  }
}
#endif /* HEMELB_COLLOIDS_LUBRICATIONBC_H */
