// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_BODYFORCES_H
#define HEMELB_COLLOIDS_BODYFORCES_H

#include "io/xml/XmlAbstractionLayer.h"
#include "units.h"
#include <map>
#include "colloids/Particle.h"

namespace hemelb
{
  namespace colloids
  {
    /** base class for all representations of a body force stored in the xml configuration file */
    class BodyForce
    {
      public:
        virtual const LatticeForceVector GetForceForParticle(const Particle&) const =0;
      protected:
        virtual ~BodyForce() {};
    };

    typedef BodyForce*(*BodyForceFactory_Create)(io::xml::XmlAbstractionLayer& xml);

    template <class TClass>
    class BodyForceFactory
    {
      public:
        static BodyForce* Create(io::xml::XmlAbstractionLayer& xml)
        {
          return TClass::ReadFromXml(xml);
        };
    };

    /** container for all body forces currently active in the simulation */
    class BodyForces
    {
      public:
        /** factory method - gets initial values from xml configuration file */
        static const void InitBodyForces(io::xml::XmlAbstractionLayer& xml);

        static const void AddBodyForce(const std::string name, const BodyForce* const);

        /** accumulates the effects of all known body forces on the particle */
        static const LatticeForceVector GetBodyForcesForParticle(const Particle& particle);

        static void ClearBodyForcesForAllSiteIds()
        {
          forceForEachSite.clear();
        }

        static void SetBodyForcesForSiteId(const site_t siteId, const LatticeForceVector force)
        {
          forceForEachSite[siteId] = force;
        }

        static const LatticeForceVector GetBodyForcesForSiteId(const site_t siteId)
        {
          return forceForEachSite[siteId];
        }

      private:
        /**
         * stores the details of all known body forces
         * the value type must be a base class pointer
         * as only pointers are type-compatible in C++
         */
        static std::map<std::string, const BodyForce* const> bodyForces;
        static std::map<site_t, LatticeForceVector> forceForEachSite;
    };
  }
}
#endif /* HEMELB_COLLOIDS_BODYFORCES_H */
