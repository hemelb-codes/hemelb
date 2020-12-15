// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_ENERGYTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_ENERGYTESTS_H

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/CellEnergy.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class EnergyTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE (EnergyTests);
          CPPUNIT_TEST (testBending);
          CPPUNIT_TEST (testVolume);
          CPPUNIT_TEST (testSurface);
          CPPUNIT_TEST (testStrain);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            BasisFixture::setUp();
            original = mesh;
            forces.resize(4, LatticeForceVector(0, 0, 0));
          }

          void tearDown()
          {
            BasisFixture::tearDown();
          }

          void testBending()
          {
            // No difference between original and current mesh
            // Hence energy is zero
            LatticeEnergy const actual0(facetBending(mesh.vertices, original, 0, 3, 1e0));
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, actual0, 1e-8);

            // Now modify mesh and check "energy" is square of angle difference
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            LatticeEnergy const actual1(facetBending(mesh.vertices, original, 0, 3, 1e0));
            mesh.vertices.back()[2] = 1e0;

            LatticeEnergy const expected(std::pow( (PI / 4e0 - std::acos(1. / std::sqrt(3.))), 2));
            CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(3.) * expected, actual1, 1e-8);
          }

          void testVolume()
          {
            // No difference between original and current mesh
            // Hence energy is zero
            LatticeEnergy const actual0(volumeEnergy(mesh.vertices, original, 1e0));
            CPPUNIT_ASSERT(helpers::is_zero(actual0));

            // Now modify mesh and check "energy" is square of volume diff
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            LatticeEnergy const actual1(volumeEnergy(mesh.vertices,
                                                     original,
                                                     2.0 * volume(original)));

            LatticeEnergy const deltaV(volume(mesh) - volume(original));
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - deltaV * deltaV));
            mesh.vertices.back()[2] = 1e0;
          }

          void testSurface()
          {
            // No difference between original and current mesh
            // Hence energy is zero
            LatticeEnergy const actual0(surfaceEnergy(mesh.vertices, original, 1e0));
            CPPUNIT_ASSERT(helpers::is_zero(actual0));

            // Now modify mesh and check "energy" is square of volume diff
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            LatticeEnergy const actual1(surfaceEnergy(mesh.vertices,
                                                      original,
                                                      2.0 * area(original)));

            LatticeEnergy const deltaS(area(mesh) - area(original));
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - deltaS * deltaS));
            mesh.vertices.back()[2] = 1e0;
          }

          void testStrain()
          {
            // No difference between original and current mesh
            // Hence energy is zero
            LatticeEnergy const actual0(strainEnergy(mesh.vertices, original, 1e0, 2e0));
            CPPUNIT_ASSERT(helpers::is_zero(actual0));

            // Now modify mesh and check "energy" is square of volume diff
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            LatticeEnergy const actual1(strainEnergy(mesh.vertices, original, 1e0, 2e0));

            LatticeEnergy const regression(0.0865562612162);
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - regression));
            mesh.vertices.back()[2] = 1e0;
          }

        protected:
          MeshData original;
          std::vector<LatticeForceVector> forces;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (EnergyTests);
    }
  }
}

#endif  // ONCE
