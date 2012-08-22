// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "colloids/BodyForces.h"
#include "colloids/BodyForceExamples.h"
#include "colloids/GraviticBodyForce.h"

namespace hemelb
{
  namespace colloids
  {
    std::map<std::string, const BodyForce* const> BodyForces::bodyForces;
    std::map<site_t, LatticeForceVector> BodyForces::forceForEachSite;

    const void BodyForces::InitBodyForces(io::xml::XmlAbstractionLayer& xml)
    {
      std::map<std::string, BodyForceFactory_Create> mapForceGenerators;
      mapForceGenerators["gravitic"] = &(GraviticBodyForceFactory::Create);
      mapForceGenerators["constant"] = &(ConstantBodyForceFactory::Create);
      mapForceGenerators["inv_r_sq"] = &(RadialBodyForceFactory::Create);

      bool ok = true;
      xml.ResetToTopLevel();
      ok &= xml.MoveToChild("colloids");
      ok &= xml.MoveToChild("bodyForces");
      if (!ok) return;

      for (std::map<std::string, BodyForceFactory_Create>::const_iterator
           iter = mapForceGenerators.begin();
           iter != mapForceGenerators.end();
           iter++)
      {
        const std::string forceClass = iter->first;
        const BodyForceFactory_Create createFunction = iter->second;
        bool found = xml.MoveToChild(forceClass);
        if (found)
        {
          while (found)
          {
            std::string forceName;
            ok &= xml.GetString("forceName", forceName);
            BodyForce* nextForce = createFunction(xml);
            BodyForces::bodyForces.insert(std::make_pair(forceName, nextForce));
            found = xml.NextSibling(forceClass);
          }
          xml.MoveToParent();
        }
      }
      xml.ResetToTopLevel();
    }

    const LatticeForceVector BodyForces::GetBodyForcesForParticle(const Particle& particle)
    {
      LatticeForceVector totalForce;
      for (std::map<std::string, const BodyForce* const>::const_iterator iter = bodyForces.begin();
           iter != bodyForces.end();
           iter++)
        totalForce += iter->second->GetForceForParticle(particle);
      return totalForce;
    }

  }
}
