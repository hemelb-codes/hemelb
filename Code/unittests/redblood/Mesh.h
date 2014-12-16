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
#include "resources/resource.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb { namespace unittests { namespace redblood {
/**
 * Class to test the simulation master.
 */
class RedBloodMeshDataIOTests : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(RedBloodMeshDataIOTests);
  CPPUNIT_TEST(testReadMesh);
  CPPUNIT_TEST(testWriteMesh);
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp() {
    std::string filename = resources::Resource("red_blood_cell.txt").Path();
    mesh = read_mesh(filename);
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

    CPPUNIT_ASSERT(any(mesh->facets.front(),  0));
    CPPUNIT_ASSERT(any(mesh->facets.front(), 12));
    CPPUNIT_ASSERT(any(mesh->facets.front(), 20));
    CPPUNIT_ASSERT(any(mesh->facets.back(), 100));
    CPPUNIT_ASSERT(any(mesh->facets.back(), 811));
    CPPUNIT_ASSERT(any(mesh->facets.back(), 244));
  }

  void testWriteMesh() {
    std::ostringstream output;
    write_mesh(output, *mesh);
    std::istringstream input(output.str());
    boost::shared_ptr<MeshData> other = read_mesh(input);
    CPPUNIT_ASSERT(other->vertices.size() == mesh->vertices.size());
    CPPUNIT_ASSERT(other->facets.size() == mesh->facets.size());
    CPPUNIT_ASSERT(compare(mesh->vertices.front() - other->vertices.front()));
    CPPUNIT_ASSERT(compare(mesh->vertices.back() - other->vertices.back()));
    for(size_t i(0); i < 3; ++i) {
      CPPUNIT_ASSERT(mesh->facets.front()[i] == other->facets.front()[i]);
      CPPUNIT_ASSERT(mesh->facets.back()[i] == other->facets.back()[i]);
    }
  }

  static bool compare(util::Vector3D<double> const &_in) {
    return _in.GetMagnitudeSquared() < 1e-8;
  }
  static bool any(boost::array<size_t, 3> const &_vec, size_t _value) {
    return _vec[0] == _value or _vec[1] == _value or _vec[2] == _value;
  }

private:
  boost::shared_ptr<MeshData> mesh;
};

class RedBloodMeshTests : public BasisFixture {

  CPPUNIT_TEST_SUITE(RedBloodMeshTests);
  CPPUNIT_TEST(testVolume);
  CPPUNIT_TEST(testBarycenter);
  CPPUNIT_TEST(testScaling);
  CPPUNIT_TEST(testTranslation);
  CPPUNIT_TEST_SUITE_END();
  public:
    void setUp() {
      BasisFixture::setUp();
      mesh.vertices.at(0) = LatticePosition(0.1, 0.1, 0.1);
    }

    void testVolume() {
      LatticePosition a(0.1, 0.1, 0.1), b(1, 0, 0), c(0, 1, 0), d(0, 0, 1);
      double const expected = std::abs((b-a).Cross(c-a).Dot(d-a)) / 6.0;
      for(size_t facet(0); facet != 4; ++facet) {
        CPPUNIT_ASSERT(std::abs(volume(mesh) - expected) < 1e-8);
        // swap vertices so orientation stays the same, check volume still
        // correct. This makes sure that volume does not depend on the order of
        // the vertices. It does depend on facet orientation. However, so do
        // the forces on the particle.
        size_t const v0(mesh.facets[facet][0]), v1(mesh.facets[facet][1]),
               v2(mesh.facets[facet][2]);
        mesh.facets[facet][0] = v1;
        mesh.facets[facet][1] = v2;
        mesh.facets[facet][2] = v0;
        CPPUNIT_ASSERT(std::abs(volume(mesh) - expected) < 1e-8);

        mesh.facets[facet][0] = v2;
        mesh.facets[facet][1] = v0;
        mesh.facets[facet][2] = v1;
        CPPUNIT_ASSERT(std::abs(volume(mesh) - expected) < 1e-8);

        // put back to first order
        mesh.facets[facet][0] = v0;
        mesh.facets[facet][1] = v1;
        mesh.facets[facet][2] = v2;
      }
    }
    void testBarycenter() {
      LatticePosition a(0.1, 0.1, 0.1), b(1, 0, 0), c(0, 1, 0), d(0, 0, 1);
      LatticePosition const expected = (a + b + c + d) * 0.25;
      CPPUNIT_ASSERT(std::abs(barycenter(mesh)[0] - expected[0]) < 1e-8);
      CPPUNIT_ASSERT(std::abs(barycenter(mesh)[1] - expected[1]) < 1e-8);
      CPPUNIT_ASSERT(std::abs(barycenter(mesh)[2] - expected[2]) < 1e-8);
    }

    void testScaling() {
      Dimensionless const scale = 2.5;
      Mesh original(mesh);
      Mesh scaled(mesh);
      scaled *= scale;

      CPPUNIT_ASSERT(
          helpers::is_zero(original.GetBarycenter() - scaled.GetBarycenter()));
      LatticePosition const first
        = (*original.GetVertices().begin() - original.GetBarycenter()) * scale
          + original.GetBarycenter();
      LatticePosition const second
        = (*(++original.GetVertices().begin()) - original.GetBarycenter()) * scale
          + original.GetBarycenter();
      CPPUNIT_ASSERT(helpers::is_zero(first - *scaled.GetVertices().begin()));
      CPPUNIT_ASSERT(
          helpers::is_zero(second - *(++scaled.GetVertices().begin())));
    }

    void testTranslation() {
      LatticePosition const offset(1, 2, 3);
      Mesh original(mesh);
      Mesh trans(mesh);
      trans += offset;

      CPPUNIT_ASSERT(helpers::is_zero(
            original.GetBarycenter() + offset - trans.GetBarycenter()
      ));
      LatticePosition const first = *original.GetVertices().begin() + offset;
      LatticePosition const second
        = *(++original.GetVertices().begin()) + offset;
      CPPUNIT_ASSERT(helpers::is_zero(first - *trans.GetVertices().begin()));
      CPPUNIT_ASSERT(
          helpers::is_zero(second - *(++trans.GetVertices().begin())));
    }
};

class TopologyTests : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(TopologyTests);
  CPPUNIT_TEST(testNodeToVertex);
  CPPUNIT_TEST(testFacetNeighbors);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    std::string filename = resources::Resource("red_blood_cube.txt").Path();
    mesh = read_mesh(filename);
    // Checks the mesh input makes sense
    CPPUNIT_ASSERT(mesh->vertices.size() == 8);
    CPPUNIT_ASSERT(mesh->facets.size() == 12);

    // Creates topology
    topo = boost::shared_ptr<MeshTopology>(new MeshTopology(*mesh));
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
      std::set<size_t> facets = topo->vertexToFacets[vertex];
      CPPUNIT_ASSERT(facets.size() == expected[vertex][0]);
      for(size_t facet(1); facet <= expected[vertex][0]; ++facet) {
        CPPUNIT_ASSERT(facets.count(expected[vertex][facet]) == 1);
      }
    }
  }

  void testFacetNeighbors() {
    CPPUNIT_ASSERT(topo->facetNeighbors.size() == 12);

    // expected[facet] = {neighbor indices}
    size_t expected[12][3] = {
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
      boost::array<size_t, 3> const &neighs
          = topo->facetNeighbors[facet];
      CPPUNIT_ASSERT(neighs.size() == 3);
      for(size_t neigh(0); neigh < 3; ++neigh)
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
  boost::shared_ptr<MeshData> mesh;
  boost::shared_ptr<MeshTopology> topo;
};

CPPUNIT_TEST_SUITE_REGISTRATION(RedBloodMeshDataIOTests);
CPPUNIT_TEST_SUITE_REGISTRATION(RedBloodMeshTests);
CPPUNIT_TEST_SUITE_REGISTRATION(TopologyTests);
}}}

#endif // ONCE
