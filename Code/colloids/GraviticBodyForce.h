// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_GRAVITICBODYFORCE_H
#define HEMELB_COLLOIDS_GRAVITICBODYFORCE_H

#include "colloids/BodyForces.h"

namespace hemelb
{
  namespace colloids
  {
    class GraviticBodyForce : public BodyForce
    {
      public:
        static BodyForce* ReadFromXml(io::xml::XmlAbstractionLayer& xml)
        {
          LatticeForceVector field;
          xml.GetDoubleVectorAndConvert("field", field);
          return new GraviticBodyForce(field);
        };

        virtual const LatticeForceVector GetForceForParticle(const Particle& particle) const
        {
          return graviticForce * particle.GetMass();
        };

      protected:
        GraviticBodyForce(const LatticeForceVector constantForce) :
          graviticForce(constantForce) {};

        const LatticeForceVector graviticForce;
    };

    class GraviticBodyForceFactory : public BodyForceFactory<GraviticBodyForce> { };
  }
}
#endif /* HEMELB_COLLOIDS_GRAVITICBODYFORCE_H */
