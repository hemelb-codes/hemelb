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
    CPPUNIT_TEST(testReadMesh);
    CPPUNIT_TEST_SUITE_END();
  public:
    void setUp() {
        std::string filename = Resource("red_blood_cell.txt").Path();
        mesh = redblood::read_mesh(filename);
    }

    void tearDown() {}

    void testReadMesh()
    {
        CPPUNIT_ASSERT(mesh);
        CPPUNIT_ASSERT(mesh->vertices.size() == 812);
        CPPUNIT_ASSERT(mesh->facets.size() == 1620);

        typedef util::Vector3D<double> Vector3d;
        Vector3d vfirst(0.850650808352040, 0.525731112119134, 0.0);
        Vector3d vlast(0.298722767882436, -0.095229509306174,
                -0.186475800096805);
        CPPUNIT_ASSERT(compare(mesh->vertices.front() - vfirst));
        CPPUNIT_ASSERT(compare(mesh->vertices.back() - vlast));

        CPPUNIT_ASSERT(mesh->facets.front().count(0) == 1);
        CPPUNIT_ASSERT(mesh->facets.front().count(12) == 1);
        CPPUNIT_ASSERT(mesh->facets.front().count(20) == 1);
        CPPUNIT_ASSERT(mesh->facets.back().count(100) == 1);
        CPPUNIT_ASSERT(mesh->facets.back().count(811) == 1);
        CPPUNIT_ASSERT(mesh->facets.back().count(244) == 1);
    }

    static bool compare(util::Vector3D<double> const &_in) {
        return _in.GetMagnitudeSquared() < 1e-8;
    }

  private:
    boost::shared_ptr<hemelb::redblood::MeshData> mesh;
};

class TopologyTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TopologyTests);
    CPPUNIT_TEST(testNodeToVertex);
    CPPUNIT_TEST(testFacetNeighbors);
    CPPUNIT_TEST_SUITE_END();

  public:
    void setUp() {
        std::string filename = Resource("red_blood_cube.txt").Path();
        mesh = redblood::read_mesh(filename);
        // Checks the mesh input makes sense
        CPPUNIT_ASSERT(mesh->vertices.size() == 8);
        CPPUNIT_ASSERT(mesh->facets.size() == 12);

        // Creates topology
        topo = boost::shared_ptr<hemelb::redblood::MeshTopology>(
            new hemelb::redblood::MeshTopology(*mesh)
        );
    }

    void tearDown() {}

    void testNodeToVertex() {
        CPPUNIT_ASSERT(topo->vertexToFacets.size() == 8);

        // expected[vertex] = {nfacets, facet indices}
        unsigned int expected[8][6] = {
            {4, 0, 1, 6, 9},
            {5, 0, 2, 3, 8, 9},
            {4, 0, 1, 3, 10},
            {5, 1, 6, 7, 10, 11},
            {5, 4, 6, 7, 8, 9},
            {4, 2, 4, 5, 8},
            {5, 2, 3, 5, 10, 11},
            {4, 4, 5, 7, 11},
        };
        for(unsigned vertex(0); vertex < 8; ++vertex) {
            std::set<unsigned int> facets = topo->vertexToFacets[vertex];
            CPPUNIT_ASSERT(facets.size() == expected[vertex][0]);
            for(unsigned facet(1); facet <= expected[vertex][0]; ++facet) {
                CPPUNIT_ASSERT(facets.count(expected[vertex][facet]) == 1);
            }
        }
    }

    void testFacetNeighbors() {
        CPPUNIT_ASSERT(topo->facetNeighbors.size() == 12);

        // expected[facet] = {nneighbors, neighbor indices}
        unsigned int expected[12][3] = {
            { 1, 3,  9},
            { 0, 6, 10},
            { 3, 5,  8},
            { 0, 2, 10},
            { 5, 7,  8},
            { 4, 2, 11},
            { 1, 7,  9},
            { 6, 4, 11},
            { 9, 4,  2},
            { 0, 6,  8},
            { 1, 3, 11},
            {10, 5,  7},
        };
        for(unsigned facet(0); facet < 12; ++facet) {
            boost::array<unsigned int, 3> const &neighs
                = topo->facetNeighbors[facet];
            CPPUNIT_ASSERT(neighs.size() == 3);
            for(unsigned neigh(0); neigh < 3; ++neigh)
                CPPUNIT_ASSERT(contains(neighs, expected[facet][neigh]));
        }
    }

    template<class T, unsigned long N>
        static bool contains(boost::array<T, N> const &_facets, T _value) {
           for(unsigned i(0); i < N; ++i)
               if(_facets[i] == _value) return true;
           return false;
        }

  private:
    boost::shared_ptr<hemelb::redblood::MeshData> mesh;
    boost::shared_ptr<hemelb::redblood::MeshTopology> topo;
};

CPPUNIT_TEST_SUITE_REGISTRATION(RedBloodMeshTests);
CPPUNIT_TEST_SUITE_REGISTRATION(TopologyTests);
}}

#endif // ONCE
