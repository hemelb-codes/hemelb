//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_DATAIO_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_DATAIO_TESTS_H

#include <sstream>
#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "resources/Resource.h"
#include "util/UnitConverter.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      /**
       * Class to test the simulation master.
       */
      class RedBloodMeshDataIOTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (RedBloodMeshDataIOTests);
          CPPUNIT_TEST (testReadMesh);
          CPPUNIT_TEST (testWriteMesh);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            std::string filename = resources::Resource("red_blood_cell.txt").Path();
            mesh = readMesh(filename);
          }

          void tearDown()
          {
          }

          void testReadMesh()
          {
            CPPUNIT_ASSERT(mesh);
            CPPUNIT_ASSERT(mesh->vertices.size() == 12);
            CPPUNIT_ASSERT(mesh->facets.size() == 20);

            typedef util::Vector3D<double> Vector3d;
            Vector3d vfirst(0.850650808352040, 0.525731112119134, 0.000000000000000);
            Vector3d vlast(0, -0.767597772408033, -0.319039022201661);
            CPPUNIT_ASSERT(compare(mesh->vertices.front() - vfirst));
            CPPUNIT_ASSERT(compare(mesh->vertices.back() - vlast));

            CPPUNIT_ASSERT(any(mesh->facets.front(), 0));
            CPPUNIT_ASSERT(any(mesh->facets.front(), 8));
            CPPUNIT_ASSERT(any(mesh->facets.front(), 4));
            CPPUNIT_ASSERT(any(mesh->facets.back(), 11));
            CPPUNIT_ASSERT(any(mesh->facets.back(), 7));
            CPPUNIT_ASSERT(any(mesh->facets.back(), 5));
          }

          void testWriteMesh()
          {
            std::ostringstream output;
            writeMesh(output, *mesh, util::UnitConverter(1, 1, LatticePosition(0, 0, 0)));
            std::istringstream input(output.str());
            std::shared_ptr<MeshData> other = readMesh(input);
            CPPUNIT_ASSERT(other->vertices.size() == mesh->vertices.size());
            CPPUNIT_ASSERT(other->facets.size() == mesh->facets.size());
            CPPUNIT_ASSERT(compare(mesh->vertices.front() - other->vertices.front()));
            CPPUNIT_ASSERT(compare(mesh->vertices.back() - other->vertices.back()));

            for (size_t i(0); i < 3; ++i)
            {
              CPPUNIT_ASSERT(mesh->facets.front()[i] == other->facets.front()[i]);
              CPPUNIT_ASSERT(mesh->facets.back()[i] == other->facets.back()[i]);
            }
          }

          static bool compare(util::Vector3D<double> const &in)
          {
            return in.GetMagnitudeSquared() < 1e-8;
          }
          static bool any(std::array<size_t, 3> const &vec, size_t value)
          {
            return vec[0] == value or vec[1] == value or vec[2] == value;
          }

        private:
          std::shared_ptr<MeshData> mesh;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (RedBloodMeshDataIOTests);
    }
  }
}

#endif  // ONCE
