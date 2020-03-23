//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_VTK_DATAIO_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_VTK_DATAIO_TESTS_H

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
       * Test reading and writing RBC meshes in VTK format(s)
       */
      class RedBloodMeshVTKDataIOTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (RedBloodMeshVTKDataIOTests);
          CPPUNIT_TEST (testVTPReadMesh);
          //CPPUNIT_TEST (testVTPWriteMesh);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            std::string filename = resources::Resource("rbc_ico_720.vtp").Path();
            mesh = readVTKMesh(filename);
          }

          void tearDown()
          {
          }

          void testVTPReadMesh()
          {
            CPPUNIT_ASSERT(mesh);
            CPPUNIT_ASSERT_EQUAL(mesh->vertices.size(), 362ul);
            CPPUNIT_ASSERT_EQUAL(mesh->facets.size(), 720ul);

            typedef util::Vector3D<double> Vector3d;
            Vector3d vfirst(0.850650808352040, 0.525731112119134, 0.000000000000000);
            Vector3d vlast(0.289241015911102, -0.180331975221633, -0.199399739503860);
            CPPUNIT_ASSERT(compare(mesh->vertices.front() - vfirst));
            CPPUNIT_ASSERT(compare(mesh->vertices.back() - vlast));

            CPPUNIT_ASSERT(any(mesh->facets.front(), 0));
            CPPUNIT_ASSERT(any(mesh->facets.front(), 12));
            CPPUNIT_ASSERT(any(mesh->facets.front(), 17));
            CPPUNIT_ASSERT(any(mesh->facets.back(), 67));
            CPPUNIT_ASSERT(any(mesh->facets.back(), 361));
            CPPUNIT_ASSERT(any(mesh->facets.back(), 157));
          }

//          void testWriteMesh()
//          {
//            std::ostringstream output;
//            writeMesh(output, *mesh, util::UnitConverter(1, 1, LatticePosition(0, 0, 0)));
//            std::istringstream input(output.str());
//            std::shared_ptr<MeshData> other = readMesh(input);
//            CPPUNIT_ASSERT(other->vertices.size() == mesh->vertices.size());
//            CPPUNIT_ASSERT(other->facets.size() == mesh->facets.size());
//            CPPUNIT_ASSERT(compare(mesh->vertices.front() - other->vertices.front()));
//            CPPUNIT_ASSERT(compare(mesh->vertices.back() - other->vertices.back()));
//
//            for (size_t i(0); i < 3; ++i)
//            {
//              CPPUNIT_ASSERT(mesh->facets.front()[i] == other->facets.front()[i]);
//              CPPUNIT_ASSERT(mesh->facets.back()[i] == other->facets.back()[i]);
//            }
//          }

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


      CPPUNIT_TEST_SUITE_REGISTRATION (RedBloodMeshVTKDataIOTests);
    }
  }
}

#endif  // ONCE
