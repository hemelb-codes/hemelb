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

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Mesh.h"

namespace hemelb { namespace unittests {
/**
 * Class to test the simulation master.
 */
class RedBloodMeshTests : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(RedBloodMeshTests);
    CPPUNIT_TEST(TestReadMesh);
    CPPUNIT_TEST_SUITE_END();
  public:
    void setUp() {
        std::string filename = Resource("red_blood_cell.txt").Path();
        data = redblood::read_mesh(filename);
    }

    void tearDown() {}

    void TestReadMesh()
    {
        CPPUNIT_ASSERT(data);
        CPPUNIT_ASSERT(data->vertices.size() == 812);
        CPPUNIT_ASSERT(data->facets.size() == 1620);

        typedef util::Vector3D<double> Vector3d;
        Vector3d vfirst(0.850650808352040, 0.525731112119134, 0.0);
        Vector3d vlast(0.298722767882436, -0.095229509306174,
                -0.186475800096805);
        CPPUNIT_ASSERT(compare(data->vertices.front() - vfirst));
        CPPUNIT_ASSERT(compare(data->vertices.back() - vlast));

        typedef util::Vector3D<unsigned int> Vector3u;
        Vector3u ffirst(0, 12, 20);
        Vector3u flast(100, 811, 244);
        CPPUNIT_ASSERT(compare(data->facets.front() - ffirst));
        CPPUNIT_ASSERT(compare(data->facets.back() - flast));
    }

    static bool compare(util::Vector3D<double> const &_in) {
        return _in.GetMagnitudeSquared() < 1e-8;
    }
    static bool compare(util::Vector3D<unsigned int> const &_in) {
        return _in.x == 0 and _in.y == 0 and _in.z == 0;
    }

  private:
    boost::shared_ptr<hemelb::redblood::MeshData> data;
};
CPPUNIT_TEST_SUITE_REGISTRATION(RedBloodMeshTests);
}}

#endif // ONCE
