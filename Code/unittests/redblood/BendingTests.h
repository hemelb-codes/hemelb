//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_BENDING_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_BENDING_TESTS_H

#include <iomanip>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "redblood/Cell.h"
#include "redblood/Mesh.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class BendingTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (BendingTests);
          CPPUNIT_TEST(testNoBendingNoNothing);
          CPPUNIT_TEST(testEnergy);
          CPPUNIT_TEST(testNumericalForces);
          CPPUNIT_TEST(testConvexities);
          CPPUNIT_TEST(testSwapFacets);
          CPPUNIT_TEST(testRotateNodeInFacet);
          CPPUNIT_TEST(testInflate);
          CPPUNIT_TEST(testMoveNodes);
          CPPUNIT_TEST(testStuff);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            using namespace hemelb::redblood;
            vertices = MeshData::Vertices{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            mesh.vertices = vertices;
            mesh.facets = MeshData::Facets{{{0, 1, 2}}, {{1, 3, 2}}};
            forces.resize(4, LatticeForceVector(0, 0, 0));
          }

          void testNoBendingNoNothing()
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, energyAndForces(), 1e-10);
            for(auto const force: forces)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(force.x, 0e0, 1e-10);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(force.y, 0e0, 1e-10);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(force.z, 0e0, 1e-10);
            }
          }
          void testEnergy()
          {
            for(auto const theta: {1e-2, 2e-2, 3e-2})
            {
              std::fill(forces.begin(), forces.end(), 0e0);
              vertices.back() = bending(mesh.vertices.back(), theta);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * theta*theta * moduli, energy(), 1e-10);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * theta*theta * moduli, energy(), 1e-10);
            }
          }

          void testStuff()
          {
            // modify template geometry to be almost planar
            auto const theta = 1e-2;
            mesh.vertices.back() = bending(LatticePosition(1, 1, 0), theta);
            // Now go to flat geometry, check energy and forces
            vertices.back() = LatticePosition(1, 1, 0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5*moduli*theta*theta, energyAndForces(), 1e-10);
            for(auto const &force: forces)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, force.x, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, force.y, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, force.z, 1e-8);
            }
          }

          void testNumericalForces()
          {
            numericalForces();
          }
          void testInflate()
          {
            for(auto &vertex: vertices)
            {
              vertex = vertex * 1.5;
            }
            testNumericalForces();
          }

          void testConvexities()
          {

            auto with_angles = [this](double theta0, double theta)
            {
              mesh.vertices.back() = bending(LatticePosition(1, 1, 0), theta0);
              vertices.back() = bending(LatticePosition(1, 1, 0), theta);
              numericalForces();
            };
            // Explores all sorts of configurations with different concavities, etc
            for(auto const theta0: {2e-1, -2e-1})
            {
              for(auto const theta: {1e-1, -1e-1, 3e-1, -3e-1})
              {
                if(std::abs(theta0 - theta) > 1e-8)
                {
                  with_angles(theta0, theta);
                }
              }
            }
          }

          void testSwapFacets()
          {
            std::swap(mesh.facets[0], mesh.facets[1]);
            testNumericalForces();
          }

          void testRotateNodeInFacet()
          {
            std::swap(mesh.facets[0][0], mesh.facets[0][1]);
            std::swap(mesh.facets[0][0], mesh.facets[0][2]);
            testNumericalForces();

            std::swap(mesh.facets[1][0], mesh.facets[1][1]);
            std::swap(mesh.facets[1][0], mesh.facets[1][2]);
            testNumericalForces();
          }

          void testMoveNodes()
          {
            for(auto &vertex: vertices)
            {
              for(int j(0); j < 6; ++j)
              {
                auto const epsilon = 1e-1;
                LatticePosition const direction(
                    j == 0 ? 1: (j == 1 ? -1: 0),
                    j == 2 ? 1: (j == 3 ? -1: 0),
                    j == 4 ? 1: (j == 5 ? -1: 0)
                );
                auto const old = vertex;
                vertex += direction * epsilon;
                numericalForces();
                vertex = old;
              }
            }
          }

          void numericalForces(LatticePosition const &direction, int node)
          {
            auto const theta0 = orientedAngle(Facet(mesh, 0), Facet(mesh, 1));
            auto const theta = orientedAngle(
                Facet(vertices, mesh.facets, 0), Facet(vertices, mesh.facets, 1));
            auto const epsilon = 1e-4;
            std::fill(forces.begin(), forces.end(), 0e0);
            auto const oldPos = vertices[node];

            auto const e0 = energyAndForces();
            if(std::abs(theta - theta0) < 1e-8)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, forces[node].Dot(direction), 1e-8);
              return;
            }

            vertices[node] += direction * epsilon;
            auto const e1 = energy();
            auto const tol= std::max(std::abs((e1 - e0)/epsilon) * 1e-2, 1e-5);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(-(e1 - e0)/epsilon, forces[node].Dot(direction), tol);
            vertices[node] = oldPos;
          };

          void numericalForces()
          {
            for(int j(0); j < 1; ++j)
            {
              LatticePosition const direction(
                  j == 0 ? 1: (j == 1 ? -1: (j >= 6 ? random() - 0.5: 0)),
                  j == 2 ? 1: (j == 3 ? -1: (j >= 6 ? random() - 0.5: 0)),
                  j == 4 ? 1: (j == 5 ? -1: (j >= 6 ? random() - 0.5: 0))
              );
              numericalForces(direction.GetNormalised(), 0);
              numericalForces(direction.GetNormalised(), 1);
              numericalForces(direction.GetNormalised(), 2);
              numericalForces(direction.GetNormalised(), 3);
            }
          }

          PhysicalEnergy energy() const
          {
            return facetBending(vertices, mesh, 0, 1, moduli);
          }
          PhysicalEnergy energyAndForces()
          {
            return facetBending(vertices, mesh, 0, 1, moduli, forces);
          }
          LatticePosition bending(LatticePosition const &v, Dimensionless theta) const
          {
            auto const &a0 = mesh.vertices[1];
            auto const &a1 = mesh.vertices[2];
            return rotationMatrix(a1 - a0, theta) * (v - a0) + a0;
          }
          LatticePosition twisting(LatticePosition const &v, Dimensionless theta) const
          {
            auto const &a0 = mesh.vertices[2];
            return rotationMatrix(LatticePosition(1, 1, 1), theta) * (v - a0) + a0;
          }

        protected:
          hemelb::redblood::MeshData mesh;
          hemelb::redblood::MeshData::Vertices vertices;
          std::vector<LatticeForceVector> forces;
          LatticeForce const moduli = 1e0;
      };


      CPPUNIT_TEST_SUITE_REGISTRATION (BendingTests);
    }
  }
}

#endif  // ONCE

