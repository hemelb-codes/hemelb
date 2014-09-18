//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARTICLEIMLP_H
#define HEMELB_UNITTESTS_REDBLOOD_PARTICLEIMLP_H

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/ParticleImpl.cc"

namespace hemelb { namespace unittests {

// Tests functionality that is *not* part of the HemeLB API
// Checks that we know how to compute geometric properties between facets
// However, HemeLB only cares about energy and forces
class FacetTests : public CppUnit::TestFixture
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
  void setUp() {
    // facets at something degrees from one another
    mesh.vertices.push_back(LatticePosition(0, 0, 0));
    mesh.vertices.push_back(LatticePosition(1, 0, 0));
    mesh.vertices.push_back(LatticePosition(0, 1, 0));
    mesh.vertices.push_back(LatticePosition(0, 0, 1));

    redblood::MeshData::t_Facet indices;
    indices[0] = 0; indices[1] = 1; indices[2] = 2;
    mesh.facets.push_back(indices);
    indices[0] = 3; indices[1] = 1; indices[2] = 2;
    mesh.facets.push_back(indices);

    main.reset(new redblood::Facet(mesh, 0));
    neighbor.reset(new redblood::Facet(mesh, 1));
  }

  void tearDown() {}

  void testNormal() {
    {
      LatticePosition const actual = normal(*main);
      LatticePosition const expected(0, 0, -1);
      CPPUNIT_ASSERT(is_zero(actual - expected));
    }
    {
      LatticePosition const actual = normal(*neighbor);
      LatticePosition const expected(-1, -1, -1);
      CPPUNIT_ASSERT(is_zero(actual - expected));
    }
  }

  void testUnitNormal() {
    {
      LatticePosition const actual = unit_normal(*main);
      LatticePosition const expected(0, 0, -1);
      CPPUNIT_ASSERT(is_zero(actual - expected));
    } {
      LatticePosition const actual = unit_normal(*neighbor);
      LatticePosition const expected(-1, -1, -1);
      CPPUNIT_ASSERT(is_zero(actual - expected.GetNormalised()));
    }
  }

  void testAngle() {
    Angle const actual0 = redblood::angle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual0 - std::acos(1.0 /std::sqrt(3.))));

    mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
    Angle const actual1 = redblood::angle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual1 - PI/4.0));
    mesh.vertices.back()[1] = 1e0;
  }

  void testCommonNodes() {
    redblood::t_IndexPair const nodes \
        = redblood::commonNodes(*main, *neighbor);
    CPPUNIT_ASSERT(nodes.first == 1 and nodes.second == 2);
    LatticePosition const edge = commonEdge(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(edge - (mesh.vertices[1] - mesh.vertices[2])));
  }
  void testSingleNodes() {
    redblood::t_IndexPair const nodes 
      = redblood::singleNodes(*main, *neighbor);
    CPPUNIT_ASSERT(nodes.first == 0 and nodes.second == 0);
  }

  void testOrientedAngle() {
    Angle const actual0 = redblood::orientedAngle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual0 - std::acos(1.0 /std::sqrt(3.))));

    // simpler angle
    mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
    Angle const actual1 = redblood::orientedAngle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual1 - PI/4.0));

    // change orientation <==> negative angle
    mesh.facets.front()[1] = 2;
    mesh.facets.front()[2] = 1;
    Angle const actual2 = redblood::orientedAngle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual2 + PI/4.0));
    mesh.facets.front()[1] = 1;
    mesh.facets.front()[2] = 2;

    mesh.vertices.back()[1] = 1e0;
  }


  static bool is_zero(util::Vector3D<double> const &_in, double _tol = 1e-8) {
    return _in.GetMagnitudeSquared() < _tol;
  }
  static bool is_zero(double const _in, double _tol = 1e-8) {
    return std::abs(_in) < _tol;
  }

protected:
  redblood::MeshData mesh;
  LatticePosition nodes[4];
  std::set<size_t> main_indices;
  std::set<size_t> neighbor_indices;
  boost::shared_ptr<redblood::Facet> main, neighbor;
};

class EnergyTests : public FacetTests {
    CPPUNIT_TEST_SUITE(EnergyTests);
    CPPUNIT_TEST(testBendingEnergy);
    CPPUNIT_TEST(testBendingForcesNoDeformation);
    CPPUNIT_TEST(testBendingForcesMoveNodes);
    CPPUNIT_TEST_SUITE_END();
  public:
    void setUp() {
      FacetTests::setUp();
      original = mesh;
      forces.resize(4, LatticeForceVector(0, 0, 0));
    }

    void tearDown() { FacetTests::tearDown(); }

    void testBendingEnergy() {
      // No difference between original and current mesh
      // Hence energy is zero
      PhysicalEnergy const actual0(
        redblood::facetBending(mesh, original, 0, 1, 1e0)
      );
      CPPUNIT_ASSERT(is_zero(actual0));

      // Now modify mesh and check "energy" is square of angle difference
      mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
      PhysicalEnergy const actual1(
        redblood::facetBending(mesh, original, 0, 1, 1e0)
      );
      mesh.vertices.back()[2] = 1e0;

      PhysicalEnergy const expected(
          std::pow((PI / 4e0 - std::acos(1./std::sqrt(3.))), 2)
      );
      CPPUNIT_ASSERT(is_zero(actual1 - expected));
    }

    void noBending(LatticePosition const &_normal, size_t _node) {
      double const epsilon(1e-5);
      mesh.vertices[_node] += _normal * epsilon;
      PhysicalEnergy const deltaE(
        redblood::facetBending(mesh, original, 0, 1, 1.0, forces)
      );
      mesh.vertices[_node] = original.vertices[_node];

      CPPUNIT_ASSERT(is_zero(deltaE / epsilon));
      for(size_t i(0); i < forces.size(); ++i)
        CPPUNIT_ASSERT(is_zero(forces[i] / epsilon));
    }

    void bending(LatticePosition const &_normal, size_t _node) {
      // Only direction perpendicular to first facet will have an effect
      double const epsilon(1e-5);
      mesh.vertices[_node] += _normal * epsilon;
      PhysicalEnergy const deltaE(
        redblood::facetBending(mesh, original, 0, 1, 1.0, forces)
      );
      mesh.vertices[_node] = original.vertices[_node];

      CPPUNIT_ASSERT(is_zero(forces[_node] + _normal * (deltaE / epsilon)));
      for(size_t i(0); i < forces.size(); ++i)
        forces[i] = LatticeForceVector(0, 0, 0);
    }

    void testBendingForcesMoveNodes() {
      // Single nodes
      bending(LatticePosition(0, 0, 1).GetNormalised(), 0);
      bending(LatticePosition(1, 1, 1).GetNormalised(), 3);
      // Common nodes
      bending(LatticePosition(0, 0, 1).GetNormalised(), 1);
      bending(LatticePosition(1, 1, 1).GetNormalised(), 1);
      bending(LatticePosition(0, 0, 1).GetNormalised(), 2);
      bending(LatticePosition(1, 1, 1).GetNormalised(), 2);

      // Single nodes
      noBending(LatticePosition(1, 0, 0).GetNormalised(), 0);
      noBending(LatticePosition(0, 1, 0).GetNormalised(), 0);
      noBending(LatticePosition(1, 1, -2).GetNormalised(), 3);
      noBending(LatticePosition(1, -1, 0).GetNormalised(), 3);
      // Common nodes
      noBending(LatticePosition(1, -1, 0).GetNormalised(), 1);
      noBending(LatticePosition(1, -1, 0).GetNormalised(), 2);
    }

    void testBendingForcesNoDeformation() {
      // No difference between original and current mesh
      // Hence energy is zero
      PhysicalEnergy const deltaE(
        redblood::facetBending(mesh, original, 0, 1, 1e0, forces)
      );
      CPPUNIT_ASSERT(is_zero(deltaE));
      for(size_t i(0); i < forces.size(); ++i)
        CPPUNIT_ASSERT(is_zero(forces[i]));
    }

  protected:
    redblood::MeshData original;
    std::vector<LatticeForceVector> forces;
};

CPPUNIT_TEST_SUITE_REGISTRATION(FacetTests);
CPPUNIT_TEST_SUITE_REGISTRATION(EnergyTests);
}}

#endif // ONCE
