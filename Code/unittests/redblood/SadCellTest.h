// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_SADCELLTEST_H
#define HEMELB_UNITTESTS_REDBLOOD_SADCELLTEST_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "redblood/Mesh.h"
#include "resources/Resource.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class SadCellTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (SadCellTests);
          CPPUNIT_TEST (test_bending);
          CPPUNIT_TEST (test_surface);
          CPPUNIT_TEST (test_volume);
          CPPUNIT_TEST (test_dilation);
          CPPUNIT_TEST (test_strain);
          CPPUNIT_TEST (test_bendingForceDirection);
          CPPUNIT_TEST (test_surfaceForceDirection);
          CPPUNIT_TEST (test_volumeForceDirection);
          CPPUNIT_TEST (test_dilationForceDirection);
          CPPUNIT_TEST (test_strainForceDirection);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            auto const normal = readMesh(resources::Resource("rbc_ico_1280.msh").Path().c_str());
            cell =
                std::make_shared<redblood::Cell>(readMesh(resources::Resource("sad.msh").Path().c_str())->vertices,
                                                 normal);
            // sad mesh is in physical units, move to something with fewer decimals.
            *cell *= 1e0 / 4.1e-6;
          }

#         define HEMELB_MACRO(name)   \
          void test_ ## name()        \
          {                           \
            cell->moduli.name = 1e0;  \
            checkNumericalForces();   \
            cell->moduli.name = 0e0;  \
          }                           \
          void test_ ## name ## ForceDirection() \
          {                                      \
            cell->moduli.name = 1e0;             \
            checkForceDirection();               \
            cell->moduli.name = 0e0;             \
          }

          HEMELB_MACRO(bending)
          ;HEMELB_MACRO(surface)
          ;HEMELB_MACRO(volume)
          ;HEMELB_MACRO(dilation)
          ;HEMELB_MACRO(strain)
          ;
#         undef HEMELB_MACRO

          void checkForceDirection()
          {
            std::vector<LatticeForceVector> forces(cell->GetVertices().size(), 0e0);
            auto const e0 = (*cell)(forces);
            for (size_t i(0); i < forces.size(); ++i)
            {
              cell->GetVertices()[i] += forces[i] * 1e-6;
            }
            auto const e1 = (*cell)();
            CPPUNIT_ASSERT(e0 > e1);
          }

          void checkNumericalForces()
          {
            std::vector<LatticeForceVector> forces(cell->GetVertices().size(), 0e0);
            (*cell)(forces);
            for (size_t i(0); i < cell->GetVertices().size(); ++i)
            {
              for (int j(0); j < 4; ++j)
              {
                auto const epsilon = 1e-3;
                LatticePosition const direction = LatticePosition(j == 0 ?
                  1 :
                  (j >= 3 ?
                    random() :
                    0),
                                                                  j == 1 ?
                                                                    1 :
                                                                    (j >= 3 ?
                                                                      random() :
                                                                      0),
                                                                  j == 2 ?
                                                                    1 :
                                                                    (j >= 3 ?
                                                                      random() :
                                                                      0)).GetNormalised();

                double const expected = derivative(direction, i, epsilon);
                double const actual = -forces[i].Dot(direction);
                double const tolerance = std::abs(actual) < 1e-7 ?
                  1e-7 :
                  1e-4 * std::abs(actual);

                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, actual, tolerance);
              }
            }
          }

          double derivative(LatticePosition const &dir, size_t i, double epsilon)
          {
            double result = 0e0;
            auto const orig = cell->GetVertices()[i];
            auto &vertex = cell->GetVertices()[i];
            std::vector<double> const coeffs { -1. / 280., 4. / 105., -1. / 5., 4. / 5. };
            for (size_t j(0); j < coeffs.size(); ++j)
            {
              vertex = orig + dir * double(coeffs.size() - j) * epsilon;
              result += coeffs[j] * (*cell)();
              vertex = orig - dir * double(coeffs.size() - j) * epsilon;
              result -= coeffs[j] * (*cell)();
            }
            vertex = orig;
            return result / epsilon;
          }

        protected:
          std::shared_ptr<redblood::Cell> cell;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (SadCellTests);
    }
  }
}

#endif  // ONCE

