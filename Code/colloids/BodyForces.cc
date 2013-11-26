// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 
#include "Exception.h"
#include "log/Logger.h"
#include "colloids/BodyForces.h"
#include "colloids/BodyForceExamples.h"
#include "colloids/GraviticBodyForce.h"

namespace hemelb
{
  namespace colloids
  {
    BodyForces* BodyForces::Load(const io::xml::Element& xml)
    {
      BodyForces* ans = new BodyForces();
      ans->Init(xml);
      return ans;
    }

    void BodyForces::Init(const io::xml::Element& bodyForcesEl)
    {
      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("Creating Body Forces.");
      if (bodyForcesEl.GetName() != "bodyForces")
        throw Exception()
            << "BodyForces::Load got XML Element with wrong name. Expected 'bodyForces', got '"
            << bodyForcesEl.GetName() << "'.";
      // Create a map from the strings used in the XML (as attribute "forceName")
      // to the factory functions that create them.

      std::map<std::string, BodyForceFactory_Create> mapForceGenerators;
      mapForceGenerators["gravitic"] = & (GraviticBodyForceFactory::Create);
      mapForceGenerators["constant"] = & (ConstantBodyForceFactory::Create);
      mapForceGenerators["inv_r_sq"] = & (RadialBodyForceFactory::Create);

      for (std::map<std::string, BodyForceFactory_Create>::const_iterator iter =
          mapForceGenerators.begin(); iter != mapForceGenerators.end(); iter++)
      {
        const std::string forceClass = iter->first;
        const BodyForceFactory_Create createFunction = iter->second;

        for (io::xml::Element forceNode = bodyForcesEl.GetChildOrNull(forceClass); forceNode
            != io::xml::Element::Missing(); forceNode = forceNode.NextSiblingOrNull(forceClass))
        {
          std::string forceName = forceNode.GetAttributeOrThrow("name");
          BodyForce* nextForce = createFunction(forceNode);
          bodyForces.insert(std::make_pair(forceName, nextForce));
        }
      }
    }

    LatticeForceVector BodyForces::GetBodyForcesForParticle(const Particle& particle) const
    {
      LatticeForceVector totalForce;
      for (std::map<std::string, const BodyForce* const >::const_iterator iter = bodyForces.begin(); iter
          != bodyForces.end(); iter++)
        totalForce += iter->second->GetForceForParticle(particle);
      return totalForce;
    }

  }
}
