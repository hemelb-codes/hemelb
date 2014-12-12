//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_GRID_AND_CELL_H
#define HEMELB_UNITTESTS_REDBLOOD_GRID_AND_CELL_H

#include <cppunit/TestFixture.h>
#include "redblood/GridAndCell.h"
#include "redblood/facet.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb { namespace unittests {

namespace {
  // Creates a mesh that is a single planar triangle
  // It still has two facets so that we can run fake forces on it.
  redblood::Mesh triangleMesh() {
    // Rotate samosa so it is in xy plane
    redblood::Mesh result = redblood::pancakeSamosa(0);
    result.GetData()->vertices[1] = LatticePosition(0, 1, 0);
    result.GetData()->vertices[2]
      = LatticePosition(std::sqrt(3.0)/2.0, 0.5, 0);
    return result;
  }
}

class CellVelocityInterpolTests : public SquareDuctTetrahedronFixture {
    CPPUNIT_TEST_SUITE(CellVelocityInterpolTests);
    CPPUNIT_TEST(testDistributionFixture);
    CPPUNIT_TEST(testLinearVelocityPerpendicularToPancakeSamosa);
    CPPUNIT_TEST(testLinearVelocityInSamosaPlane);
    CPPUNIT_TEST_SUITE_END();

    typedef lb::lattices::D3Q15 D3Q15;
    typedef lb::kernels::LBGK<D3Q15> Kernel;
public:

    void setUp() { SquareDuctTetrahedronFixture::setUp(); }

    // Checks fixture function do what they should do
    void testDistributionFixture();

    // Linear velociy profile on grid, perpendicular to the samosa
    // The interpolated velocities should be constant across the samosa
    void testLinearVelocityPerpendicularToPancakeSamosa();
    // Linear velociy profile on grid, in-plane with the samosa
    // The interpolated velocities should evolve linearly with respect to the
    // input gradient.
    void testLinearVelocityInSamosaPlane();

protected:
    // Creates a mesh that is a single planar triangle
    // It still has two facets so that we can run fake forces on it.
    virtual redblood::Mesh initial_mesh() const {
      return redblood::pancakeSamosa(0);
    }
    virtual size_t refinement() const { return 3; }
};


class CellForceSpreadTests : public SquareDuctTetrahedronFixture {
    CPPUNIT_TEST_SUITE(CellForceSpreadTests);
      CPPUNIT_TEST(testIsZeroFarFromMembrane);
      CPPUNIT_TEST(testIsSymmetric);
      CPPUNIT_TEST(testIsIncreasing);
      CPPUNIT_TEST(testIsLinear);
    CPPUNIT_TEST_SUITE_END();

    typedef lb::lattices::D3Q15 D3Q15;
    typedef lb::kernels::LBGK<D3Q15> Kernel;
public:

    void setUp();
    void setUpForces();

    // Moves fixture so barycenter is at given position and computes spread
    // force at origin
    LatticeForceVector force_at_center(LatticePosition const &_position);
    void testIsZeroFarFromMembrane();
    void testIsSymmetric();
    void testIsIncreasing();
    void testIsLinear();

  protected:
    std::vector<LatticeForceVector> forces;
    LatticePosition direction, intensity, center;
    virtual redblood::Mesh initial_mesh() const { return triangleMesh(); }
    virtual size_t refinement() const { return 3; }
};

class CellForceSpreadWithWallTests : public SquareDuctTetrahedronFixture {
    CPPUNIT_TEST_SUITE(CellForceSpreadWithWallTests);
      CPPUNIT_TEST(testAPIAssumptions);
      CPPUNIT_TEST(testNode2WallCutoff);
    CPPUNIT_TEST_SUITE_END();
    typedef lb::lattices::D3Q15 D3Q15;
    typedef lb::kernels::LBGK<D3Q15> Kernel;
  public:
    void testAPIAssumptions();
    void testNode2WallCutoff();
  protected:
    int siteID(LatticeVector const &) const;
    LatticeVector GetSolidWall() const;
    virtual redblood::Mesh initial_mesh() const { return triangleMesh(); }
    virtual size_t refinement() const { return 3; }
};


// Sets up a linear velocity profile
// Uses a macro so we can define a number of variables in one go
# define HEMELB_LINEAR_VELOCITY_PROFILE(AX, AY, AZ)                        \
  LatticePosition const gradient= LatticePosition(AX, AY, AZ).Normalise(); \
  /* any number big enough to avoid negative populations */                \
  Dimensionless const non_negative_population(CubeSize() * 3);             \
  helpers::Linear linear(non_negative_population, gradient);               \
  helpers::Linear linear_inv(2 * non_negative_population, -gradient);      \
  helpers::setUpDistribution<D3Q15>(latDat, 0, linear);                    \
  helpers::setUpDistribution<D3Q15>(latDat, 1, linear_inv)


void CellVelocityInterpolTests::testDistributionFixture() {
  helpers::ZeroOutFOld(latDat);

  HEMELB_LINEAR_VELOCITY_PROFILE(2., 4., 6.);
  // Test assumes static pop at index == 0 as assumed by macro
  CPPUNIT_ASSERT_DOUBLES_EQUAL(D3Q15::CX[0], 0e0, 1e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(D3Q15::CY[0], 0e0, 1e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(D3Q15::CZ[0], 0e0, 1e-8);

  size_t const N = 4;
  LatticeVector const a[N] = {
     LatticeVector(2, 4, 3),
     LatticeVector(10, 16, 5),
     LatticeVector(20, 3, 10),
     LatticeVector(22, 8, 15)
  };
  for(size_t i(0); i < N; ++i) {
    size_t const index = latDat->GetContiguousSiteId(a[i]);
    LatticePosition const pos(a[i][0], a[i][1], a[i][2]);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(
        linear(pos),
        latDat->GetSite(index).GetFOld<D3Q15>()[0],
        1e-8
    );
    CPPUNIT_ASSERT_DOUBLES_EQUAL(
        linear_inv(pos),
        latDat->GetSite(index).GetFOld<D3Q15>()[1],
        1e-8
    );
    for(size_t j(2); j < D3Q15::NUMVECTORS; ++j) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(
          0,
          latDat->GetSite(index).GetFOld<D3Q15>()[j],
          1e-8
      );
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(
        3.0 * non_negative_population,
        latDat->GetSite(index).GetFOld<D3Q15>()[0]
        + latDat->GetSite(index).GetFOld<D3Q15>()[1],
        1e-8
    );
  }
}

void CellVelocityInterpolTests :: testLinearVelocityPerpendicularToPancakeSamosa() {
  // direction perpendicular to plane
  helpers::ZeroOutFOld(latDat);
  LatticePosition const normal(redblood::Facet(*mesh.GetData(), 0).normal());
  HEMELB_LINEAR_VELOCITY_PROFILE(normal.x, normal.y, normal.z);

  // Perform interpolation
  std::vector<LatticePosition> displacements;
  redblood::velocitiesOnMesh<Kernel>(
      mesh, *latDat, redblood::stencil::FOUR_POINT, displacements);

  // Compute expected velocities
  typedef std::vector<LatticePosition> :: const_iterator const_iterator;
  const_iterator i_disp = displacements.begin();
  const_iterator const i_end = displacements.end();
  LatticePosition const expected(*i_disp);
  for(++i_disp; i_disp != i_end; ++i_disp) {
    CPPUNIT_ASSERT(helpers::is_zero(*i_disp - expected));
  }
}

void CellVelocityInterpolTests :: testLinearVelocityInSamosaPlane() {
  // Figures out an in-plane direction
  helpers::ZeroOutFOld(latDat);
  redblood::Facet const shapeFacet(*mesh.GetData(), 0);
  LatticePosition const inplane(shapeFacet.edge(0) + shapeFacet.edge(1) * 0.5);
  HEMELB_LINEAR_VELOCITY_PROFILE(inplane.x, inplane.y, inplane.z);

  // Perform interpolation
  std::vector<LatticePosition> displacements;
  redblood::velocitiesOnMesh<Kernel>(
      mesh, *latDat, redblood::stencil::FOUR_POINT, displacements);

  // Computes what the interpolation should be
  typedef std::vector<LatticePosition> :: const_iterator const_iterator;
  LatticeDistance const
    x0 = gradient.Dot(mesh.GetVertex(0)),
    x1 = gradient.Dot(mesh.GetVertex(1));
  PhysicalVelocity const
    v0 = displacements[0],
    v1 = displacements[1];
  redblood::MeshData::t_Vertices::const_iterator
    i_vertex(mesh.GetVertices().begin() + 2);
  const_iterator i_disp = displacements.begin() + 2;
  const_iterator const i_end = displacements.end();
  for(; i_disp != i_end; ++i_disp, ++i_vertex) {
    PhysicalVelocity const expected(
        (v0 - v1) * ((i_vertex->Dot(gradient) - x1) / (x0 - x1)) + v1
    );
    CPPUNIT_ASSERT(helpers::is_zero(*i_disp - expected));
  }
}
# undef HEMELB_LINEAR_VELOCITY_PROFILE

LatticeForceVector CellForceSpreadTests :: force_at_center(
    LatticePosition const &_position) {
  mesh += _position - mesh.GetBarycenter();
  helpers::ZeroOutForces(latDat);
  redblood::details::spreadForce2Grid(
      mesh, redblood::details::SpreadForces(forces, *latDat),
      redblood::stencil::FOUR_POINT
  );
  return latDat->GetSite(center).GetForce();
}


void CellForceSpreadTests :: setUp() {
  center = LatticePosition(CubeSize() / 2, CubeSize() / 2, CubeSize() / 2);
  SquareDuctTetrahedronFixture::setUp();
  setUpForces();
}
void CellForceSpreadTests :: setUpForces() {
  direction = LatticePosition(1, 2, 3);
  intensity = LatticePosition(3, 2, 1).Normalise();
  forces.resize(mesh.GetNumberOfNodes());
  typedef redblood::MeshData::t_Vertices::const_iterator const_iterator;
  const_iterator i_vertex = mesh.GetVertices().begin();
  const_iterator const i_end = mesh.GetVertices().end();
  std::vector<LatticeForceVector> :: iterator i_force = forces.begin();
  for(; i_vertex != i_end; ++i_vertex, ++i_force)
    *i_force = direction * i_vertex->Dot(intensity);
}

void CellForceSpreadTests :: testIsZeroFarFromMembrane() {
  // Very far away from pancake samosa
  LatticeForceVector const faraway = force_at_center(
      center + LatticePosition(0, 0, 3)
  );
  CPPUNIT_ASSERT(helpers::is_zero(faraway));

  LatticeForceVector const justToFar = force_at_center(
      center + LatticePosition(0, 0, 2)
  );
  CPPUNIT_ASSERT(helpers::is_zero(justToFar));

  LatticeForceVector const justInside = force_at_center(
      center + LatticePosition(0, 0, 2.0 - 1e-4)
  );
  CPPUNIT_ASSERT(not helpers::is_zero(justInside));
  CPPUNIT_ASSERT(justInside[0] > 1e-8 and justInside[0] < 1e-4);
}


void CellForceSpreadTests :: testIsSymmetric() {
  size_t const N(10);
  for(size_t i(0); i < N; ++i) {
    LatticePosition const displacement(
      0, 0, 2e0 - Dimensionless(i * 2) / Dimensionless(N)
    );
    LatticeForceVector const left = force_at_center(center + displacement);
    LatticeForceVector const right = force_at_center(center + displacement);
    CPPUNIT_ASSERT(helpers::is_zero(left - right));
  }
}

void CellForceSpreadTests :: testIsIncreasing() {
  size_t const N(10);
  LatticeForceVector last(0, 0, 0);
  for(size_t i(0); i < N; ++i) {
    LatticePosition const displacement(
      0, 0, 2e0 - Dimensionless(i * 2 + 1) / Dimensionless(N)
    );
    LatticeForceVector const current = force_at_center(center + displacement);
    CPPUNIT_ASSERT(current[0] > last[0]);
    CPPUNIT_ASSERT(current[1] > last[1]);
    CPPUNIT_ASSERT(current[2] > last[2]);
    last = current;
  }
}

void CellForceSpreadTests :: testIsLinear() {
  size_t const N(5);
  mesh = redblood::Cell(redblood::refine(mesh, 4));
  setUpForces();
  // x0, x1 should be further than 2 from the edges
  // Only linear if samosa appears as infinite plane
  // with sufficiently dense vertices
  LatticePosition const
    x0(center[0], center[1] - 0.5, center[2] - 0.1),
    x1(center[0], center[1] + 0.5, center[2] - 0.1);
  LatticeForceVector const
    v0(force_at_center(x0)),
    v1(force_at_center(x1));

  LatticeForceVector const a(
      (v1 - v0) / (direction.Dot(x1) - direction.Dot(x0))
  );
  Dimensionless const tolerance(
      std::max(
        std::max(std::abs((v0 - v1)[0]), std::abs((v0 - v1)[1])),
        std::abs((v0 - v1)[2])
      ) * 1e-3
  );

  for(size_t i(0); i < N; ++i) {
    LatticePosition const x
        = (x1 - x0) * (Dimensionless(i + 1) / Dimensionless(N + 2)) + x0;
    LatticeForceVector const expected(
        a * (direction.Dot(x) - direction.Dot(x0)) + v0
    );
    CPPUNIT_ASSERT(helpers::is_zero(expected - force_at_center(x), tolerance));
  }
}

// IsValidLatticeSite doesn't always return false when site is invalid
// It seems to be more of a quick and innacurate check
int CellForceSpreadWithWallTests :: siteID(
    LatticeVector const &_position) const {
  site_t siteid;
  proc_t procid;
  bool const valid = latDat->GetContiguousSiteId(_position, procid, siteid);
  return valid ? siteid: -1;
}
LatticeVector CellForceSpreadWithWallTests :: GetSolidWall() const {
  return latDat->GetGlobalSiteMaxes()
     - LatticeVector(CubeSize() / 2, -1, CubeSize() / 2);
}

void CellForceSpreadWithWallTests :: testAPIAssumptions() {
  // Nearest solid node. Does not have associated data.
  LatticeVector const solid = GetSolidWall();
  // Nearest wet wall node. Has fluid particles and is closest node to wall.
  LatticeVector const wetwall = solid + LatticeVector(0, -1, 0);
  LatticeVector const left = wetwall + LatticeVector(-1, 0, 0);
  LatticeVector const right = wetwall + LatticeVector(1, 0, 0);
  LatticeVector const fluid = wetwall + LatticeVector(0, -1, 0);

  // Checking I know where walls and fluid sites are.
  CPPUNIT_ASSERT(siteID(solid) < 0);
  CPPUNIT_ASSERT(siteID(wetwall) >= 0);
  CPPUNIT_ASSERT(siteID(left) >= 0);
  CPPUNIT_ASSERT(siteID(right) >= 0);
  CPPUNIT_ASSERT(not latDat->GetSite(siteID(fluid)).IsSolid());
  CPPUNIT_ASSERT(not latDat->GetSite(siteID(wetwall)).IsSolid());
  CPPUNIT_ASSERT(not latDat->GetSite(siteID(left)).IsSolid());
  CPPUNIT_ASSERT(not latDat->GetSite(siteID(right)).IsSolid());
  CPPUNIT_ASSERT(not latDat->GetSite(siteID(fluid)).IsWall());
  CPPUNIT_ASSERT(latDat->GetSite(siteID(wetwall)).IsWall());
  CPPUNIT_ASSERT(latDat->GetSite(siteID(left)).IsWall());
  CPPUNIT_ASSERT(latDat->GetSite(siteID(right)).IsWall());
}

void CellForceSpreadWithWallTests :: testNode2WallCutoff() {
  // Fluid sites next to wall
  LatticeVector const solid = GetSolidWall();
  LatticeVector const wetwall = solid + LatticeVector(0, -1, 0);
  LatticeVector const left = wetwall + LatticeVector(-1, 0, 0);
  LatticeVector const right = wetwall + LatticeVector(1, 0, 0);

  // Makes sure there are no force except node-Wall interaction
  helpers::ZeroOutForces(latDat);
  mesh.moduli = redblood::Cell::Moduli();
  mesh.nodeWall.intensity = 1.;
  // Min distance to wall, so we can check cutoff
  helpers::SetWallDistance(latDat, 0.3);

  // Each test case comes with cutoff distance and mesh position
  PhysicalDistance const cutoffs[] = {
    0.25,
    0.35,
    0.35,
    std::sqrt(0.2 * 0.2 + 2) + 0.1,
    -1, // Stops loop!!!
  };
  LatticePosition const positions[] = {
    wetwall.cast<PhysicalDistance>(),
    wetwall.cast<PhysicalDistance>(),
    wetwall.cast<PhysicalDistance>() - LatticePosition(0, 0.2, 0),
    wetwall.cast<PhysicalDistance>() - LatticePosition(0, 0.2, 0)
  };
  bool const expected[] = {
    true, true, true,
    false, true, true,
    true, true, true,
    false, false, false,
  };

  // Loop over test cases breaks on special marker (negative) cutoff
  for(size_t i(0); cutoffs[i] > 0e0; ++i) {
    helpers::ZeroOutForces(latDat);
    mesh += positions[i] - mesh.GetVertex(0);
    mesh.nodeWall.cutoff = cutoffs[i];

    redblood::forcesOnGridWithWallInteraction<D3Q15>(
        mesh,  *latDat, redblood::stencil::FOUR_POINT);

    bool const atWall = helpers::is_zero(latDat->GetSite(wetwall).GetForce());
    bool const atLeft = helpers::is_zero(latDat->GetSite(left).GetForce());
    bool const atRight = helpers::is_zero(latDat->GetSite(right).GetForce());
    CPPUNIT_ASSERT(atWall == expected[3*i]);
    CPPUNIT_ASSERT(atLeft == expected[3*i + 1]);
    CPPUNIT_ASSERT(atRight == expected[3*i + 2]);
  }
}

CPPUNIT_TEST_SUITE_REGISTRATION(CellVelocityInterpolTests);
CPPUNIT_TEST_SUITE_REGISTRATION(CellForceSpreadTests);
CPPUNIT_TEST_SUITE_REGISTRATION(CellForceSpreadWithWallTests);
}}

#endif // ONCE

