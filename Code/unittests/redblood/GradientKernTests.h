// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_GRADIENTKERNTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_GRADIENTKERNTESTS_H

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      // Tests certain node movement that do not result in forces/energy
      class GradientKernTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE (GradientKernTests);
          CPPUNIT_TEST (testBending);
          CPPUNIT_TEST (testVolume);
          CPPUNIT_TEST (testSurface);CPPUNIT_TEST_SUITE_END();

#define HEMELB_MACRO(NAME, FUNCTION_CALL)                                   \
  void NAME(LatticePosition const &normal, size_t node)                     \
  {                                                                         \
    std::vector<LatticeForceVector> forces(4, LatticeForceVector(0, 0, 0)); \
    double const epsilon(1e-5);                                             \
    MeshData newmesh(mesh);                                                 \
    newmesh.vertices[node] += normal * epsilon;                             \
    LatticeEnergy const deltaE(FUNCTION_CALL);                              \
    CPPUNIT_ASSERT(helpers::is_zero(deltaE / epsilon));                     \
    for (size_t i(0); i < forces.size(); ++i)                               \
      CPPUNIT_ASSERT(helpers::is_zero(forces[i]));                          \
  }

          HEMELB_MACRO(noBending, facetBending(newmesh.vertices, mesh, 0, 3, 1.0, forces))
          ;HEMELB_MACRO(noVolume, volumeEnergy(newmesh.vertices, mesh, 1.0, forces))
          ;HEMELB_MACRO(noSurface, surfaceEnergy(newmesh.vertices, mesh, 1.0, forces))
          ;
#undef HEMELB_MACRO

          void testBending()
          {
            // Single nodes
            noBending(LatticeForceVector(1, 0, 0).GetNormalised(), 0);
            noBending(LatticeForceVector(0, 1, 0).GetNormalised(), 0);
            noBending(LatticeForceVector(1, 1, -2).GetNormalised(), 3);
            noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 3);
            // Common nodes
            noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 1);
            noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 2);
          }

          void testVolume()
          {
            noVolume(LatticePosition(-1, 1, 0).GetNormalised(), 0);
            noVolume(LatticePosition(1, 1, -2).GetNormalised(), 0);
            noVolume(LatticePosition(0, 1, 0).GetNormalised(), 1);
            noVolume(LatticePosition(0, 0, 1).GetNormalised(), 1);
            noVolume(LatticePosition(1, 0, 0).GetNormalised(), 2);
            noVolume(LatticePosition(0, 0, 1).GetNormalised(), 2);
            noVolume(LatticePosition(1, 0, 0).GetNormalised(), 3);
            noVolume(LatticePosition(0, 1, 0).GetNormalised(), 3);
          }

          void testSurface()
          {
            // Get rid of some surfaces first, to simplify kernel
            MeshData const save(mesh);
            mesh.facets.resize(2);
            mesh.facets.back() = save.facets.back();

            try
            {
              noSurface(LatticePosition(0, 0, 1).GetNormalised(), 0);
              noSurface(LatticePosition(1, 1, 1).GetNormalised(), 3);
              noSurface(LatticePosition(-1, 1, 0).GetNormalised(), 0);
              noSurface(LatticePosition(-1, 1, 0).GetNormalised(), 3);
            }
            catch (...)
            {
              mesh = save;
              throw;
            }

            mesh = save;
          }
      };


      CPPUNIT_TEST_SUITE_REGISTRATION (GradientKernTests);
    }
  }
}

#endif  // ONCE
