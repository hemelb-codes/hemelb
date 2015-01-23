//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_GRADIENT_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_GRADIENT_TESTS_H

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      struct EnergyFunctional
      {
          EnergyFunctional(MeshData const &original) :
              original(original)
          {
          }
          ;
          MeshData const &original;
      };
#define HEMELB_OP_MACRO(CLASS, FUNCTION)                                                           \
  struct CLASS : public EnergyFunctional                                                           \
  {                                                                                                \
    CLASS(MeshData const &original) : EnergyFunctional(original)                                   \
    {                                                                                              \
    }                                                                                              \
    PhysicalEnergy operator()(MeshData const &mesh) const                                          \
    {                                                                                              \
      return FUNCTION(mesh.vertices, original, 1.0);                                               \
    }                                                                                              \
    PhysicalEnergy operator()(MeshData const &mesh, std::vector<LatticeForceVector> &forces) const \
    {                                                                                              \
      return FUNCTION(mesh.vertices, original, 1.0, forces);                                       \
    }                                                                                              \
  }
      HEMELB_OP_MACRO(VolumeEnergy, volumeEnergy);
      HEMELB_OP_MACRO(SurfaceEnergy, surfaceEnergy);
#undef HEMELB_OP_MACRO

      struct BendingEnergy : public EnergyFunctional
      {
          BendingEnergy(MeshData const &original) :
              EnergyFunctional(original)
          {
          }
          PhysicalEnergy operator()(MeshData const &mesh) const
          {
            return facetBending(mesh.vertices, original, 0, 1, 1.0)
                + facetBending(mesh.vertices, original, 0, 2, 1.0)
                + facetBending(mesh.vertices, original, 0, 3, 1.0)
                + facetBending(mesh.vertices, original, 1, 2, 1.0)
                + facetBending(mesh.vertices, original, 1, 3, 1.0)
                + facetBending(mesh.vertices, original, 2, 3, 1.0);
          }
          PhysicalEnergy operator()(MeshData const &mesh,
                                    std::vector<LatticeForceVector> &forces) const
          {
            return facetBending(mesh.vertices, original, 0, 1, 1.0, forces)
                + facetBending(mesh.vertices, original, 0, 2, 1.0, forces)
                + facetBending(mesh.vertices, original, 0, 3, 1.0, forces)
                + facetBending(mesh.vertices, original, 1, 2, 1.0, forces)
                + facetBending(mesh.vertices, original, 1, 3, 1.0, forces)
                + facetBending(mesh.vertices, original, 2, 3, 1.0, forces);
          }
      };

      struct StrainEnergy : public EnergyFunctional
      {
          StrainEnergy(MeshData const &original) :
              EnergyFunctional(original)
          {
          }
          PhysicalEnergy operator()(MeshData const &mesh) const
          {
            return strainEnergy(mesh.vertices, original, 1.0, 2.0);
          }
          PhysicalEnergy operator()(MeshData const &mesh,
                                    std::vector<LatticeForceVector> &forces) const
          {
            return strainEnergy(mesh.vertices, original, 1.0, 2.0, forces);
          }
      };

      class GradientTests : public EnergyVsGradientFixture
      {
          CPPUNIT_TEST_SUITE (GradientTests);
          CPPUNIT_TEST (testBending);
          CPPUNIT_TEST (testSurface);
          CPPUNIT_TEST (testVolume);
          CPPUNIT_TEST (testStrain);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            BasisFixture::setUp();
            original = mesh;
            mesh.vertices[0] += LatticePosition(-0.01, 0.02342, 0.03564);
            mesh.vertices[1] += LatticePosition(0.0837, -0.012632, 0.0872935);
            mesh.vertices[2] += LatticePosition(0.02631, -0.00824223, -0.098362);
          }

          void testBending()
          {
            energyVsForces(BendingEnergy(original));
          }
          void testSurface()
          {
            energyVsForces(SurfaceEnergy(original));
          }
          void testVolume()
          {
            energyVsForces(VolumeEnergy(original));
          }
          void testStrain()
          {
            energyVsForces(StrainEnergy(original));
          }

        protected:
          MeshData original;
      };


      CPPUNIT_TEST_SUITE_REGISTRATION (GradientTests);
    }
  }
}

#endif  // ONCE
