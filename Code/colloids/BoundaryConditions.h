// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_COLLOIDS_BOUNDARYCONDITIONS_H
#define HEMELB_COLLOIDS_BOUNDARYCONDITIONS_H

#include "units.h"
#include <vector>
#include <map>
#include "colloids/Particle.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace io
  {
    namespace xml
    {
      class Element;
    }
  }

  namespace colloids
  {
    /** base class for all representations of a boundary condition */
    class BoundaryCondition
    {
      public:
        /** applies the boundary condition by directly modifying the particle
         *  returns false if the particle must now be deleted, true otherwise
         */
        virtual const bool DoSomethingToParticle(Particle&, const std::vector<LatticePosition>) =0;
        virtual const std::vector<Particle> CreateNewParticles() =0;
      protected:
        virtual ~BoundaryCondition()
        {
        }
        ;
    };

    typedef BoundaryCondition*(*BoundaryConditionFactory_Create)(io::xml::Element& xml);

    template<class TClass>
    class BoundaryConditionFactory
    {
      public:
        static BoundaryCondition* Create(io::xml::Element& xml)
        {
          return TClass::ReadFromXml(xml);
        }
        ;
    };

    /** container for all colloids boundary conditions currently active in the simulation */
    class BoundaryConditions
    {
      public:
        /** factory method - gets initial values from xml configuration file */
        static BoundaryConditions* Load(const geometry::LatticeData& latDat,
                                                          const io::xml::Element& bcEl);

        void AddBoundaryCondition(const std::string name, const BoundaryCondition* const );

        /** allows all registered boundary conditions to do something to this particle
         returns false if the particle should be deleted, true otherwise
         */
        bool DoSomeThingsToParticle(const LatticeTimeStep currentTimestep,
                                    Particle& particle);

        /** allows all registered boundary conditions to create new particles */
        const std::vector<Particle> CreateNewParticles();

      private:
        BoundaryConditions(const geometry::LatticeData& latDat);
        /**
         * stores the details of all known body forces
         * the value type must be a base class pointer
         * as only pointers are type-compatible in C++
         */
        std::vector<BoundaryCondition*> boundaryConditionsWall;
        std::vector<BoundaryCondition*> boundaryConditionsInlet;
        std::vector<BoundaryCondition*> boundaryConditionsOutlet;

        const geometry::LatticeData& latticeData;
    };
  }
}
#endif /* HEMELB_COLLOIDS_BOUNDARYCONDITIONS_H */
