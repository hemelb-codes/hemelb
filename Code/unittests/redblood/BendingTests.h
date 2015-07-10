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
          CPPUNIT_TEST(testTwisting);
          CPPUNIT_TEST(testLargeAngles);
          CPPUNIT_TEST(testCommonEdge);
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

          void testTwisting()
          {
            for(auto const theta: {0e0, 1e-2, 2e-2, 3e-2})
            {
              for(auto const phi: {1e-2, 2e-2, 3e-2})
              {
                std::fill(forces.begin(), forces.end(), 0e0);
                vertices.back() = bending(twisting(mesh.vertices.back(), phi), theta);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * theta*theta * moduli, energy(), 1e-10);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * theta*theta * moduli, energy(), 1e-10);
              }
            }
          }

          void testLargeAngles()
          {
            // modify template geometry to be almost planar
            auto const theta = 1e-2;
            mesh.vertices.back() = bending(LatticePosition(1, 1, 0), theta);

            // check we still get zero if vertices == template
            vertices.back() = mesh.vertices.back();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, energyAndForces(), 1e-10);
            for(auto const force: forces)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(force.x, 0e0, 1e-10);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(force.y, 0e0, 1e-10);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(force.z, 0e0, 1e-10);
            }

            // Now go to flat geometry, check energy and forces
            vertices.back() = LatticePosition(1, 1, 0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5*moduli*theta*theta, energyAndForces(), 1e-10);
            CPPUNIT_ASSERT(forces[0].z < -1e-4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(forces[0].z, forces[3].z, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(
                forces[0].z, forces[0].Dot(LatticePosition(0, 0, 1)), 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(
                forces[3].z, forces[3].Dot(LatticePosition(0, 0, 1)), 1e-8);

            // Then go to other concavity and check result
            vertices.back() = bending(LatticePosition(1, 1, 0), -theta);
            std::fill(forces.begin(), forces.end(), 0e0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0*moduli*theta*theta, energyAndForces(), 1e-10);
            CPPUNIT_ASSERT(forces[0].z < -1e-4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(
                forces[0].z, forces[0].Dot(LatticePosition(0, 0, 1)), 1e-8);
            CPPUNIT_ASSERT(forces[3].z < -1e-4);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(forces[0].GetMagnitude(), forces[3].GetMagnitude(), 1e-8);
          }

          // Checks direction of restoration force on common edge
          void testCommonEdge()
          {
            auto const dir1 = (LatticePosition(0, 0, 0.5) - vertices[1]).GetNormalised();
            auto const dir2 = (LatticePosition(0, 0, 0.5) - vertices[2]).GetNormalised();
            vertices[1] += dir1 * 0.01;
            vertices[2] += dir2 * 0.01;
            energyAndForces();
            CPPUNIT_ASSERT(forces[1].Dot(dir1) < 0e0);
            CPPUNIT_ASSERT(forces[2].Dot(dir2) < 0e0);
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

