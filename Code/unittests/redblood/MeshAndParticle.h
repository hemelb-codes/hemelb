//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_MESH_AND_PARTICLE_H
#define HEMELB_UNITTESTS_REDBLOOD_MESH_AND_PARTICLE_H

#include <cppunit/TestFixture.h>
#include "redblood/MeshAndParticle.h"
#include "redblood/facet.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb { namespace unittests {


class ParticleVelocityInterpolTests : public SquareDuctTetrahedronFixture {
    CPPUNIT_TEST_SUITE(ParticleVelocityInterpolTests);
    CPPUNIT_TEST(testDistributionFixture);
    CPPUNIT_TEST(testLinearVelocityPerpendicularToPancakeSamosa);
    CPPUNIT_TEST(testLinearVelocityInSamosaPlane);
    CPPUNIT_TEST_SUITE_END();

    typedef lb::lattices::D3Q15 D3Q15;
    typedef lb::kernels::LBGK<D3Q15> Kernel;
public:

    void setUp() {
      SquareDuctTetrahedronFixture::setUp();
    }

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


void ParticleVelocityInterpolTests::testDistributionFixture() {
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

void ParticleVelocityInterpolTests :: testLinearVelocityPerpendicularToPancakeSamosa() {
  // direction perpendicular to plane
  helpers::ZeroOutFOld(latDat);
  LatticePosition const normal(redblood::Facet(*mesh.GetData(), 0).normal());
  HEMELB_LINEAR_VELOCITY_PROFILE(normal.x, normal.y, normal.z);

  // Perform interpolation
  std::vector<LatticePosition> displacements;
  redblood::compute_displacement<Kernel>(
      mesh, *latDat, redblood::stencil::FOUR_POINT, displacements);

  // Compute expected velocities
  typedef std::vector<LatticePosition> :: const_iterator const_iterator;
  const_iterator i_disp = displacements.begin();
  const_iterator const i_end = displacements.end();
  LatticePosition const expected(*i_disp);
  for(++i_disp; i_disp != i_end; ++i_disp) {
    CPPUNIT_ASSERT(is_zero(*i_disp - expected));
  }
}

void ParticleVelocityInterpolTests :: testLinearVelocityInSamosaPlane() {
  // Figures out an in-plane direction
  helpers::ZeroOutFOld(latDat);
  redblood::Facet const shapeFacet(*mesh.GetData(), 0);
  LatticePosition const inplane(shapeFacet.edge(0) + shapeFacet.edge(1) * 0.5);
  HEMELB_LINEAR_VELOCITY_PROFILE(inplane.x, inplane.y, inplane.z);

  // Perform interpolation
  std::vector<LatticePosition> displacements;
  redblood::compute_displacement<Kernel>(
      mesh, *latDat, redblood::stencil::FOUR_POINT, displacements);

  // Computes what the interpolation should be
  typedef std::vector<LatticePosition> :: const_iterator const_iterator;
  LatticeDistance const
    x0 = gradient.Dot(mesh.GetVertices()[0]),
    x1 = gradient.Dot(mesh.GetVertices()[1]);
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
    CPPUNIT_ASSERT(is_zero(*i_disp - expected));
  }
}
# undef HEMELB_LINEAR_VELOCITY_PROFILE

CPPUNIT_TEST_SUITE_REGISTRATION(ParticleVelocityInterpolTests);
}}

#endif // ONCE

