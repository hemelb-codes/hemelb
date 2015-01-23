//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARTICLEIMPL_H
#define HEMELB_UNITTESTS_REDBLOOD_PARTICLEIMPL_H

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Cell.impl.cc"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {

      // Tests functionality that is *not* part of the HemeLB API
      // Checks that we know how to compute geometric properties between facets
      // However, HemeLB only cares about energy and forces
      class FacetTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE(FacetTests);
          CPPUNIT_TEST(testNormal);
          CPPUNIT_TEST(testUnitNormal);
          CPPUNIT_TEST(testAngle);
          CPPUNIT_TEST(testCommonNodes);
          CPPUNIT_TEST(testOrientedAngle);
          CPPUNIT_TEST(testSingleNodes);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            BasisFixture::setUp();
            main.reset(new Facet(mesh, 0));
            neighbor.reset(new Facet(mesh, 3));
          }

          void tearDown() {}

          void testNormal()
          {
            {
              LatticePosition const actual = main->normal();
              LatticePosition const expected(0, 0, -1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected));
            }
            {
              LatticePosition const actual = neighbor->normal();
              LatticePosition const expected(1, 1, 1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected));
            }
          }

          void testUnitNormal()
          {
            {
              LatticePosition const actual = main->unitNormal();
              LatticePosition const expected(0, 0, -1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected));
            }
            {
              LatticePosition const actual = neighbor->unitNormal();
              LatticePosition const expected(1, 1, 1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected.GetNormalised()));
            }
          }

          void testAngle()
          {
            Angle const actual0 = angle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual0 - std::acos(-1.0 / std::sqrt(3.))));

            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            Angle const actual1 = angle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - 3.0 * PI / 4.0));
            mesh.vertices.back()[1] = 1e0;
          }

          void testCommonNodes()
          {
            IndexPair const nodes = commonNodes(*main, *neighbor);
            CPPUNIT_ASSERT(nodes.first == 1 and nodes.second == 2);
            LatticePosition const edge = commonEdge(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(
                             edge - (mesh.vertices[1] - mesh.vertices[2])
                           ));
          }
          void testSingleNodes()
          {
            IndexPair const nodes = singleNodes(*main, *neighbor);
            // indices into the list of vertices of each facet,
            // not into list of all vertices.
            CPPUNIT_ASSERT(nodes.first == 0 and nodes.second == 1);
          }

          void testOrientedAngle()
          {
            Angle const actual0 = orientedAngle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual0 + std::acos(-1.0 / std::sqrt(3.))));

            // simpler angle
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            Angle const actual1 = orientedAngle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual1 + 3.0 * PI / 4.0));

            // change orientation <==> negative angle
            mesh.facets.front()[1] = 2;
            mesh.facets.front()[2] = 1;
            Angle const actual2 = orientedAngle(Facet(mesh, 0), Facet(mesh, 3));
            CPPUNIT_ASSERT(helpers::is_zero(actual2 + PI / 4.0));
            mesh.facets.front()[1] = 1;
            mesh.facets.front()[2] = 2;

            mesh.vertices.back()[1] = 1e0;
          }
        protected:
          LatticePosition nodes[4];
          std::set<size_t> main_indices;
          std::set<size_t> neighbor_indices;
          std::shared_ptr<Facet> main, neighbor;
      };

      class EnergyTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE(EnergyTests);
          CPPUNIT_TEST(testBending);
          CPPUNIT_TEST(testVolume);
          CPPUNIT_TEST(testSurface);
          CPPUNIT_TEST(testStrain);
          CPPUNIT_TEST_SUITE_END();
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
            PhysicalEnergy const actual0(
              facetBending(mesh.vertices, original, 0, 3, 1e0)
            );
            CPPUNIT_ASSERT(helpers::is_zero(actual0));

            // Now modify mesh and check "energy" is square of angle difference
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            PhysicalEnergy const actual1(
              facetBending(mesh.vertices, original, 0, 3, 1e0)
            );
            mesh.vertices.back()[2] = 1e0;

            PhysicalEnergy const expected(
              std::pow((PI / 4e0 - std::acos(1. / std::sqrt(3.))), 2)
            );
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - expected));
          }

          void testVolume()
          {
            // No difference between original and current mesh
            // Hence energy is zero
            PhysicalEnergy const actual0(volumeEnergy(mesh.vertices, original, 1e0));
            CPPUNIT_ASSERT(helpers::is_zero(actual0));

            // Now modify mesh and check "energy" is square of volume diff
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            PhysicalEnergy const actual1(
              volumeEnergy(mesh.vertices, original, 2.0 * volume(original))
            );

            PhysicalEnergy const deltaV(volume(mesh) - volume(original));
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - deltaV * deltaV));
            mesh.vertices.back()[2] = 1e0;
          }

          void testSurface()
          {
            // No difference between original and current mesh
            // Hence energy is zero
            PhysicalEnergy const actual0(
              surfaceEnergy(mesh.vertices, original, 1e0)
            );
            CPPUNIT_ASSERT(helpers::is_zero(actual0));

            // Now modify mesh and check "energy" is square of volume diff
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            PhysicalEnergy const actual1(
              surfaceEnergy(mesh.vertices, original, 2.0 * surface(original))
            );

            PhysicalEnergy const deltaS(surface(mesh) - surface(original));
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - deltaS * deltaS));
            mesh.vertices.back()[2] = 1e0;
          }

          void testStrain()
          {
            // No difference between original and current mesh
            // Hence energy is zero
            PhysicalEnergy const actual0(
              strainEnergy(mesh.vertices, original, 1e0, 2e0)
            );
            CPPUNIT_ASSERT(helpers::is_zero(actual0));

            // Now modify mesh and check "energy" is square of volume diff
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            PhysicalEnergy const actual1(
              strainEnergy(mesh.vertices, original, 1e0, 2e0)
            );

            PhysicalEnergy const regression(0.0865562612162);
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - regression));
            mesh.vertices.back()[2] = 1e0;
          }

        protected:
          MeshData original;
          std::vector<LatticeForceVector> forces;
      };



      struct EnergyFunctional
      {
        EnergyFunctional(MeshData const &original) : original(original) {};
        MeshData const &original;
      };
#define HEMELB_OP_MACRO(CLASS, FUNCTION)                       \
  struct CLASS : public EnergyFunctional {                     \
    CLASS(MeshData const &original)                           \
      : EnergyFunctional(original) {}                         \
    PhysicalEnergy operator()(MeshData const &mesh) const {   \
      return FUNCTION(mesh.vertices, original, 1.0);          \
    }                                                          \
    PhysicalEnergy operator()(MeshData const &mesh,           \
                              std::vector<LatticeForceVector> &forces) const{         \
      return FUNCTION(mesh.vertices, original, 1.0, forces); \
    }                                                          \
  }
      HEMELB_OP_MACRO(VolumeEnergy, volumeEnergy);
      HEMELB_OP_MACRO(SurfaceEnergy, surfaceEnergy);
#undef HEMELB_OP_MACRO

      struct BendingEnergy : public EnergyFunctional
      {
        BendingEnergy(MeshData const &original)
          : EnergyFunctional(original) {}
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
        StrainEnergy(MeshData const &original) : EnergyFunctional(original) {}
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
          CPPUNIT_TEST_SUITE(GradientTests);
          CPPUNIT_TEST(testBending);
          CPPUNIT_TEST(testSurface);
          CPPUNIT_TEST(testVolume);
          CPPUNIT_TEST(testStrain);
          CPPUNIT_TEST_SUITE_END();
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

      // Tests certain node movement that do not result in forces/energy
      class GradientKernTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE(GradientKernTests);
          CPPUNIT_TEST(testBending);
          CPPUNIT_TEST(testVolume);
          CPPUNIT_TEST(testSurface);
          CPPUNIT_TEST_SUITE_END();

#   define HEMELB_MACRO(NAME, FUNCTION_CALL)                                  \
  void NAME(LatticePosition const &normal, size_t node) {                 \
    std::vector<LatticeForceVector> forces(4, LatticeForceVector(0, 0, 0)); \
    double const epsilon(1e-5);                                             \
    MeshData newmesh(mesh);                                                 \
    newmesh.vertices[node] += normal * epsilon;                           \
    PhysicalEnergy const deltaE(FUNCTION_CALL);                             \
    CPPUNIT_ASSERT(helpers::is_zero(deltaE / epsilon));                     \
    for(size_t i(0); i < forces.size(); ++i)                                \
      CPPUNIT_ASSERT(helpers::is_zero(forces[i]));                          \
  }

          HEMELB_MACRO(noBending,
                       facetBending(newmesh.vertices, mesh, 0, 3, 1.0, forces));
          HEMELB_MACRO(noVolume,
                       volumeEnergy(newmesh.vertices, mesh, 1.0, forces));
          HEMELB_MACRO(noSurface,
                       surfaceEnergy(newmesh.vertices, mesh, 1.0, forces));
#   undef HEMELB_MACRO

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
            catch(...)
            {
              mesh = save;
              throw;
            }

            mesh = save;
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(FacetTests);
      CPPUNIT_TEST_SUITE_REGISTRATION(EnergyTests);
      CPPUNIT_TEST_SUITE_REGISTRATION(GradientTests);
      CPPUNIT_TEST_SUITE_REGISTRATION(GradientKernTests);
    }
  }
}

#endif // ONCE
