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
     *  However, the calculation should take into account the buoyancy caused by the fluid
     *  so the "mass" of the particle should really be the difference between its density
     *  and the density of the surrounding fluid, which can be calculated from viscosity?
     *  This needs changes to PersistedParticle, config.xml and MPI types (mass->density)
     */
    class GraviticBodyForce : public BodyForce
    {
      public:
        static BodyForce* ReadFromXml(const io::xml::Element& xml)
        {
          LatticeAccelerationVector g;
          // TODO: convert to lattice units
          io::xml::Element gEl = xml.GetChildOrThrow("acceleration");
          configuration::GetDimensionalValue(gEl, "m s^-2", g);

          return new GraviticBodyForce(g);
        };

        virtual const LatticeForceVector GetForceForParticle(const Particle& particle) const
        {
          return accelerationDueToGravity * particle.GetMass();
        };
        const LatticeAccelerationVector& GetAccelerationDueToGravity() const
        {
          return accelerationDueToGravity;
        }
      protected:
        GraviticBodyForce(const LatticeForceVector constantForce) :
          accelerationDueToGravity(constantForce) {};

        const LatticeAccelerationVector accelerationDueToGravity;
    };

    class GraviticBodyForceFactory : public BodyForceFactory<GraviticBodyForce> { };
  }
}
#endif /* HEMELB_COLLOIDS_GRAVITICBODYFORCE_H */
