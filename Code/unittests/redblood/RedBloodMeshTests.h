//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_MESH_H
#define HEMELB_UNITTESTS_REDBLOOD_MESH_H

#include <sstream>
#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "resources/Resource.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      // Tests geometry properties
      class RedBloodMeshTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE (RedBloodMeshTests);
          CPPUNIT_TEST (testOrientation);
          CPPUNIT_TEST (testStandardTestMeshOrientation);
          CPPUNIT_TEST (testVolume);
          CPPUNIT_TEST (testBarycenter);
          CPPUNIT_TEST (testScaling);
          CPPUNIT_TEST (testTranslation);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            BasisFixture::setUp();
            mesh.vertices.at(0) = LatticePosition(0.1, 0.1, 0.1);
          }

          void checkMeshOrientation(MeshData const &mesh)
          {
            MeshData copy(mesh);
            orientFacets(copy);
            for(size_t i(0); i < mesh.facets.size(); ++i)
            {
              CPPUNIT_ASSERT_EQUAL(mesh.facets[i][0], copy.facets[i][0]);
              CPPUNIT_ASSERT_EQUAL(mesh.facets[i][1], copy.facets[i][1]);
              CPPUNIT_ASSERT_EQUAL(mesh.facets[i][2], copy.facets[i][2]);
            }
          }

          void testOrientation()
          {
            checkMeshOrientation(mesh);

            // change facet order of copy for next test
            MeshData copy(mesh);
            for(auto &facet: copy.facets)
            {
              std::swap(facet[0], facet[2]);
            }
            // re-orient
            orientFacets(copy);
            checkMeshOrientation(copy);
          }

          void testStandardTestMeshOrientation()
          {
            checkMeshOrientation(*readMesh(Resource("red_blood_cell.txt").Path()));
            checkMeshOrientation(*tetrahedron(2).GetData());
            checkMeshOrientation(*icoSphere(5).GetData());
          }

          void testVolume()
          {
            LatticePosition a(0.1, 0.1, 0.1), b(1, 0, 0), c(0, 1, 0), d(0, 0, 1);
            double const expected = std::abs( (b - a).Cross(c - a).Dot(d - a)) / 6.0;

            for (size_t facet(0); facet != 4; ++facet)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, volume(mesh), 1e-8);
              // swap vertices so orientation stays the same, check volume still
              // correct. This makes sure that volume does not depend on the order of
              // the vertices. It does depend on facet orientation. However, so do
              // the forces on the particle.
              size_t const v0(mesh.facets[facet][0]), v1(mesh.facets[facet][1]),
                  v2(mesh.facets[facet][2]);
              mesh.facets[facet][0] = v1;
              mesh.facets[facet][1] = v2;
              mesh.facets[facet][2] = v0;
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, volume(mesh), 1e-8);

              mesh.facets[facet][0] = v2;
              mesh.facets[facet][1] = v0;
              mesh.facets[facet][2] = v1;
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, volume(mesh), 1e-8);

              // put back to first order
              mesh.facets[facet][0] = v0;
              mesh.facets[facet][1] = v1;
              mesh.facets[facet][2] = v2;
            }
          }
          void testBarycenter()
          {
            LatticePosition a(0.1, 0.1, 0.1), b(1, 0, 0), c(0, 1, 0), d(0, 0, 1);
            LatticePosition const expected = (a + b + c + d) * 0.25;
            CPPUNIT_ASSERT(std::abs(barycenter(mesh)[0] - expected[0]) < 1e-8);
            CPPUNIT_ASSERT(std::abs(barycenter(mesh)[1] - expected[1]) < 1e-8);
            CPPUNIT_ASSERT(std::abs(barycenter(mesh)[2] - expected[2]) < 1e-8);
          }

          void testScaling()
          {
            Dimensionless const scale = 2.5;
            Mesh original(mesh);
            Mesh scaled(mesh);
            scaled *= scale;

            CPPUNIT_ASSERT(helpers::is_zero(original.GetBarycenter() - scaled.GetBarycenter()));
            LatticePosition const first = (*original.GetVertices().begin()
                - original.GetBarycenter()) * scale + original.GetBarycenter();
            LatticePosition const second = (* (++original.GetVertices().begin())
                - original.GetBarycenter()) * scale + original.GetBarycenter();
            CPPUNIT_ASSERT(helpers::is_zero(first - *scaled.GetVertices().begin()));
            CPPUNIT_ASSERT(helpers::is_zero(second - * (++scaled.GetVertices().begin())));
          }

          void testTranslation()
          {
            LatticePosition const offset(1, 2, 3);
            Mesh original(mesh);
            Mesh trans(mesh);
            trans += offset;

            CPPUNIT_ASSERT(helpers::is_zero(original.GetBarycenter() + offset
                - trans.GetBarycenter()));
            LatticePosition const first = *original.GetVertices().begin() + offset;
            LatticePosition const second = * (++original.GetVertices().begin()) + offset;
            CPPUNIT_ASSERT(helpers::is_zero(first - *trans.GetVertices().begin()));
            CPPUNIT_ASSERT(helpers::is_zero(second - * (++trans.GetVertices().begin())));
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (RedBloodMeshTests);
    }
  }
}

#endif  // ONCE
