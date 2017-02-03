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
    /** THIS CLASS SHOULD NOT BE USED IN ITS CURRENT FORM
     *
     *  It is intended to approximate the external force of gravity on colloid particles.
     *  However, the calculation should take into account the boyancy caused by the fluid
     *  so the "mass" of the particle should really be the difference between its density
     *  and the density of the surrounding fluid, which can be calculated from viscosity?
     *  This needs changes to PersistedParticle, config.xml and MPI types (mass->density)
     */
    class GraviticBodyForce : public BodyForce
    {
      public:
        static BodyForce* ReadFromXml(io::xml::Element& xml)
        {
          LatticeForceVector field;
          // TODO: convert to lattice units
          io::xml::Element fieldElem = xml.GetChildOrThrow("field");
          fieldElem.GetAttributeOrThrow("x", field.x);
          fieldElem.GetAttributeOrThrow("y", field.y);
          fieldElem.GetAttributeOrThrow("z", field.z);

          return new GraviticBodyForce(field);
        }
        ;

        virtual const LatticeForceVector GetForceForParticle(const Particle& particle) const
        {
          return graviticForce * particle.GetMass();
        }
        ;

      protected:
        GraviticBodyForce(const LatticeForceVector constantForce) :
            graviticForce(constantForce)
        {
        }
        ;

        const LatticeForceVector graviticForce;
    };

    class GraviticBodyForceFactory : public BodyForceFactory<GraviticBodyForce>
    {
    };
  }
}
#endif /* HEMELB_COLLOIDS_GRAVITICBODYFORCE_H */
