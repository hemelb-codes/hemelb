// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_REDBLOODMESHVTKDATAIOTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_REDBLOODMESHVTKDATAIOTESTS_H

#include <sstream>
#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "resources/Resource.h"
#include "util/UnitConverter.h"
#include "unittests/redblood/Fixtures.h"
#include <vtkPolyData.h>

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
          CPPUNIT_TEST (testOrientOriginalRBCMesh);
          CPPUNIT_TEST (testOrientTimmSimRBCMesh);
          CPPUNIT_TEST_SUITE_END();

	  redblood::VTKMeshIO io = {};

        public:
          void setUp()
          {
          }

          void tearDown()
          {
          }

          void testVTPReadMesh()
          {
            std::string filename = resources::Resource("rbc_ico_720.vtp").Path();
            std::shared_ptr<MeshData> mesh = io.readFile(filename, true);

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

          void testOrientOriginalRBCMesh()
          {
            // This file has 684 out of 720 faces oriented inwards
            std::string filename = resources::Resource("rbc_ico_720.vtp").Path();

            std::shared_ptr<MeshData> meshData;
            vtkSmartPointer<vtkPolyData> polyData;
            std::tie(meshData, polyData) = io.readUnoriented(redblood::VTKMeshIO::Mode::file, filename);

            auto numSwaps = orientFacets(*meshData, *polyData);

            CPPUNIT_ASSERT_EQUAL(numSwaps, 684u);
            CPPUNIT_ASSERT(volume(*meshData) > 0.0);
          }

          void testOrientTimmSimRBCMesh()
          {
            // This file has all 720 faces oriented inwards
            std::string filename = resources::Resource("992Particles_rank3_26_t992.vtp").Path();

            std::shared_ptr<MeshData> meshData;
            vtkSmartPointer<vtkPolyData> polyData;
            std::tie(meshData, polyData) = io.readUnoriented(redblood::VTKMeshIO::Mode::file, filename);

            auto numSwaps = orientFacets(*meshData, *polyData);

            CPPUNIT_ASSERT_EQUAL(numSwaps, 0u);
            CPPUNIT_ASSERT(volume(*meshData) > 0.0);
          }

          static bool compare(util::Vector3D<double> const &in)
          {
            return in.GetMagnitudeSquared() < 1e-8;
          }
          static bool any(std::array<redblood::IdType, 3> const &vec, size_t value)
          {
            return vec[0] == value or vec[1] == value or vec[2] == value;
          }

      };


      CPPUNIT_TEST_SUITE_REGISTRATION (RedBloodMeshVTKDataIOTests);
    }
  }
}

#endif  // ONCE
