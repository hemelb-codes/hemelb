//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_INTERPOLATION_H
#define HEMELB_UNITTESTS_REDBLOOD_INTERPOLATION_H

#include <cppunit/TestFixture.h>
#include "redblood/interpolation.h"
#include "redblood/VelocityInterpolation.h"
#include "lb/lattices/D3Q15.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/Comparisons.h"

namespace hemelb { namespace unittests {

class InterpolationTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(InterpolationTests);
    CPPUNIT_TEST(testIndexIterator);
    CPPUNIT_TEST(testOffLattice);
    CPPUNIT_TEST(testOffLatticeZeroOutsideStencil);
    CPPUNIT_TEST(testInterpolateLinearFunction);
    CPPUNIT_TEST(testInterpolateQuadraticFunction);
    CPPUNIT_TEST_SUITE_END();

    struct PlanarFunction {
      LatticePosition operator()(LatticePosition const &_pos) const {
        return LatticePosition(
          LatticePosition(1, 1, 1).Dot(_pos),
          LatticePosition(-1, 2, 1).Dot(_pos),
          LatticePosition(0, 0, 1).Dot(_pos)
        );
      }
      LatticePosition operator()(Dimensionless _x, Dimensionless _y,
          Dimensionless _z) const {
        return operator()(LatticePosition(_x, _y, _z));
      }
    };

    struct QuadraticFunction {
      LatticePosition operator()(LatticePosition const &_pos) const {
        Dimensionless const offset(0);
        return LatticePosition(
          (LatticePosition(1, 1, 1).Dot(_pos) - offset)
          * (LatticePosition(1, 1, 1).Dot(_pos) - offset),
          (LatticePosition(0, 1, 0).Dot(_pos) - offset)
          * (LatticePosition(0, 1, 0).Dot(_pos) - offset),
          (LatticePosition(0, 0, 1).Dot(_pos) - offset)
          * (LatticePosition(0, 0, 1).Dot(_pos) - offset)
        );
      }
      LatticePosition operator()(Dimensionless _x, Dimensionless _y,
          Dimensionless _z) const {
        return operator()(LatticePosition(_x, _y, _z));
      }
    };

public:

    void testIndexIterator() {
      using hemelb::redblood::IndexIterator;

      LatticeVector vectors[] = {
        LatticeVector(4, 3, 2),
        LatticeVector(4, 3, 3),
        LatticeVector(4, 3, 4),
        LatticeVector(4, 4, 2),
        LatticeVector(4, 5, 2),
        LatticeVector(5, 3, 2),
        LatticeVector(6, 3, 2),
        LatticeVector(6, 5, 4),
      };
      size_t incs[] = {
        0, 1, 1, 1, 3, 3, 9, 8,
        666 // break
      };

      // Checks iteration goes through correct sequence
      IndexIterator iterator(LatticeVector(5, 4, 3), 1);
      for(size_t i(0); incs[i] < 666; ++i)  {
        for(size_t j(0); j < incs[i]; ++j, ++iterator);
        CPPUNIT_ASSERT(helpers::is_zero(*iterator - vectors[i]));
        CPPUNIT_ASSERT(iterator.isValid());
      }
      // Checks iterator becomes invalid
      ++iterator;
      CPPUNIT_ASSERT(not iterator.isValid());
    }

    void testOffLattice() {
      using hemelb::redblood::InterpolationIterator;

      LatticePosition const pos(56.51, 52.9, 15.2);
      hemelb::redblood::stencil::FourPoint stencil;
      InterpolationIterator iterator(pos, stencil);

      LatticeVector vectors[] = {
        LatticeVector(55, 51, 14),
        LatticeVector(55, 51, 15),
        LatticeVector(55, 51, 16),
        LatticeVector(55, 51, 17),
        LatticeVector(55, 52, 14),
        LatticeVector(56, 51, 14),
        LatticeVector(58, 54, 17),
      };
      size_t incs[] = {
        0, 1, 1, 1, 1, 12, 47,
        666 // break
      };

      // Checks iteration goes through correct sequence
      for(size_t i(0); incs[i] < 666; ++i)  {
        for(size_t j(0); j < incs[i]; ++j, ++iterator);
        CPPUNIT_ASSERT(helpers::is_zero(*iterator - vectors[i]));
        CPPUNIT_ASSERT(helpers::is_zero(
              stencil(pos - vectors[i]) - iterator.weight()
        ));
        CPPUNIT_ASSERT(iterator.isValid());
      }
      ++iterator;
      CPPUNIT_ASSERT(not iterator.isValid());
    }

    void testOffLatticeZeroOutsideStencil() {
      using hemelb::redblood::InterpolationIterator;
      LatticePosition const pos(56.51, 52.9, 15.2);
      hemelb::redblood::stencil::FourPoint stencil;
      InterpolationIterator iterator(pos, stencil);
      // Checks that outside iteration box, weights are zero
      LatticeVector zero_vecs[] = {
        LatticeVector(57, 53, 13),
        LatticeVector(57, 53, 18),
        LatticeVector(57, 50, 17),
        LatticeVector(57, 55, 17),
        LatticeVector(54, 53, 17),
        LatticeVector(59, 53, 17)
      };
      for(size_t i(0); i < 6; ++i) {
        LatticeVector const dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
        CPPUNIT_ASSERT(helpers::is_zero(stencil(pos - zero_vecs[i])));
        // checks we are one step outside the iteration box only.
        // this is really a test on the zero_vecs data, eg a test of the test.
        size_t const one_non_zero =
            size_t(not helpers::is_zero(stencil(pos + dx - zero_vecs[i])))
          + size_t(not helpers::is_zero(stencil(pos - dx - zero_vecs[i])))
          + size_t(not helpers::is_zero(stencil(pos + dy - zero_vecs[i])))
          + size_t(not helpers::is_zero(stencil(pos - dy - zero_vecs[i])))
          + size_t(not helpers::is_zero(stencil(pos + dz - zero_vecs[i])))
          + size_t(not helpers::is_zero(stencil(pos - dz - zero_vecs[i])));
        CPPUNIT_ASSERT(one_non_zero == 1);
      }
    }

    template<class FUNCTION> void check(Dimensionless _x, Dimensionless _y,
        Dimensionless _z, Dimensionless _tolerance = 1e-8) {
      using hemelb::redblood::interpolate;
      using hemelb::redblood::stencil::FOUR_POINT;

      FUNCTION func;
      LatticePosition expected(func(_x, _y, _z));
      LatticePosition actual(interpolate(func, _x, _y, _z, FOUR_POINT));
      CPPUNIT_ASSERT(helpers::is_zero(actual - expected, _tolerance));
    }

    // Test interpolation when the point is on the grid
    void testInterpolateLinearFunction() {
      check<PlanarFunction>(0, 0, 0);
      check<PlanarFunction>(0.1, 0.5, 0.6);
      check<PlanarFunction>(-5.1, 0.5, 8.7);
      check<PlanarFunction>(-5, 0, -1);
    }

    void testInterpolateQuadraticFunction() {
      QuadraticFunction quad;
      // Error depends on variation on scale larger than stencil
      Dimensionless const tolerance(
          (quad(0, 0, 0) - quad(12, 12, 12)).GetMagnitude() * 1e-2);
      check<QuadraticFunction>(0, 0, 0, tolerance);
      check<QuadraticFunction>(0.1, 0.5, 0.6, tolerance);
      check<QuadraticFunction>(-5.1, 0.5, 8.7, tolerance);
      check<QuadraticFunction>(-5, 0, -1, tolerance);
    }
};

class VelocityInterpolationTests :
  public helpers::FourCubeBasedTestFixture {
    typedef lb::lattices::D3Q15 t_Lattice;
    typedef lb::kernels::LBGK<t_Lattice> t_Kernel;
    CPPUNIT_TEST_SUITE(VelocityInterpolationTests);
    CPPUNIT_TEST(testLatticeDataFunctor);
    CPPUNIT_TEST_SUITE_END();

  public:
    void setUp() {
      FourCubeBasedTestFixture::setUp();

      LatticeVector const min(latDat->GetGlobalSiteMins());
      LatticeVector const max(latDat->GetGlobalSiteMaxes());
      for(size_t i(min[0]); i <= max[0]; ++i)
        for(size_t j(min[1]); j <= max[1]; ++j)
          for(size_t k(min[2]); k <= max[2]; ++k) {
            distribn_t fold[t_Lattice::NUMVECTORS];
            for(size_t u(0); u < t_Lattice::NUMVECTORS; ++u) fold[u] = 0e0;
            fold[2] = i; fold[4] = j; fold[6] = k;
            site_t const local(
                latDat->GetContiguousSiteId(LatticeVector(i, j, k))
            );
            latDat->SetFOld<t_Lattice>(local, fold);
          }
    }
    void tearDown() { FourCubeBasedTestFixture::tearDown(); }

    void testLatticeDataFunctor() {
      LatticeVector sites[] = {
        LatticeVector(1, 1, 1),
        LatticeVector(2, 1, 1),
        LatticeVector(1, 2, 1),
        LatticeVector(1, 1, 2),
        LatticeVector(4, 4, 4),
        LatticeVector(0, 0, 0)
      };

      typedef lb::lattices::D3Q15 t_Lattice;
      typedef lb::kernels::LBGK<t_Lattice> t_Kernel;
      redblood::details::VelocityFromLatticeData<t_Kernel> const functor(*latDat);
      for(size_t i(0); sites[i][0] != 0; ++i) {
        LatticeVelocity const expected =
          LatticeVelocity(
              t_Lattice::CX[2], t_Lattice::CY[2], t_Lattice::CZ[2]
          ) * Dimensionless(sites[i][0])
          + LatticeVelocity(
              t_Lattice::CX[4], t_Lattice::CY[4], t_Lattice::CZ[4]
          ) * Dimensionless(sites[i][1])
          + LatticeVelocity(
              t_Lattice::CX[6], t_Lattice::CY[6], t_Lattice::CZ[6]
          ) * Dimensionless(sites[i][2]);
        LatticeVelocity const actual = functor(sites[i]);
        CPPUNIT_ASSERT(helpers::is_zero(actual - expected));
      }
    }
};


CPPUNIT_TEST_SUITE_REGISTRATION(VelocityInterpolationTests);
}}

#endif // ONCE

