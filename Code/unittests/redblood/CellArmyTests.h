//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLARMY_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLARMY_H

#include <cppunit/TestFixture.h>
#include "unittests/redblood/Fixtures.h"
#include "redblood/CellArmy.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      //! Mock cell for ease of use
      class FakeCell : public hemelb::redblood::Cell
      {
        public:
          mutable size_t nbcalls = 0;
#         ifndef CPP11_HAS_CONSTRUCTOR_INHERITANCE
	      FakeCell(Mesh const &mesh)
                  : Cell(mesh)
              {
              }
	      FakeCell(std::shared_ptr<MeshData> const &mesh)
                  : Cell(mesh)
              {
              }
#         else
              using hemelb::redblood::Cell::Cell;
#         endif 
          //! Facet bending energy
          virtual PhysicalEnergy operator()() const override
          {
            return 0;
          }
          //! Facet bending energy
          virtual PhysicalEnergy operator()(std::vector<LatticeForceVector> &) const override
          {
            ++nbcalls;
            return 0;
          }
          virtual LatticeForceVector WallInteractionForce(
              LatticePosition const &vertex, LatticePosition const &wall) const
          {
            return 0;
          }
      };

      class CellArmyTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE (CellArmyTests);
          CPPUNIT_TEST (testCell2Fluid);
          CPPUNIT_TEST (testCell2FluidWithoutCells);
          CPPUNIT_TEST (testCellInsertion);
          CPPUNIT_TEST (testCellRemoval);
          CPPUNIT_TEST (testCellOutput);
          CPPUNIT_TEST (testFluid2Cell);CPPUNIT_TEST_SUITE_END();

          PhysicalDistance const cutoff = 5.0;
          PhysicalDistance const halo = 2.0;
          typedef lb::lattices::D3Q15 D3Q15;
          typedef lb::kernels::GuoForcingLBGK<lb::lattices::D3Q15> Kernel;

        public:
          void testCell2Fluid();
          void testCell2FluidWithoutCells();
          void testFluid2Cell();
          void testCellInsertion();
          void testCellRemoval();
          void testCellOutput();

        private:
          virtual size_t CubeSize() const
          {
            return 32 + 2;
          }
      };

      void CellArmyTests::testCell2FluidWithoutCells()
      {
        CellContainer cells;
        redblood::CellArmy<Kernel> army(*latDat, cells, cutoff, halo);
        army.cell2Cell.cutoff = 0.5;
        army.cell2Cell.intensity = 1.0;
        army.Cell2FluidInteractions();
      }

      void CellArmyTests::testCell2Fluid()
      {
        // Fixture all pairs far from one another
        auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
        assert(cells.size() == 2);
        assert(std::dynamic_pointer_cast<FakeCell>((*cells.begin()))->nbcalls == 0);
        assert(std::dynamic_pointer_cast<FakeCell>((*std::next(cells.begin())))->nbcalls == 0);

        helpers::ZeroOutFOld(latDat);
        helpers::ZeroOutForces(latDat);

        redblood::CellArmy<Kernel> army(*latDat, cells, cutoff, halo);
        army.cell2Cell.cutoff = 0.5;
        army.cell2Cell.intensity = 1.0;
        army.Cell2FluidInteractions();

        CPPUNIT_ASSERT(std::dynamic_pointer_cast<FakeCell>((*cells.begin()))->nbcalls == 1);
        CPPUNIT_ASSERT(std::dynamic_pointer_cast<FakeCell>((*std::next(cells.begin())))->nbcalls == 1);

        for (size_t i(0); i < latDat->GetLocalFluidSiteCount(); ++i)
        {
          CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(i).GetForce()));
        }

        LatticePosition const n0(15 - 0.1, 15.5, 15.5);
        LatticePosition const n1(15 + 0.2, 15.5, 15.5);
        (*cells.begin())->GetVertices().front() = n0;
        (*std::next(cells.begin()))->GetVertices().front() = n1;
        army.updateDNC();
        army.Cell2FluidInteractions();

        CPPUNIT_ASSERT(std::dynamic_pointer_cast<FakeCell>((*cells.begin()))->nbcalls == 2);
        CPPUNIT_ASSERT(std::dynamic_pointer_cast<FakeCell>((*std::next(cells.begin())))->nbcalls == 2);
        CPPUNIT_ASSERT(not helpers::is_zero(latDat->GetSite(15, 15, 15).GetForce()));
      }

      void CellArmyTests::testCellInsertion()
      {
        auto cell = std::make_shared<FakeCell>(pancakeSamosa());

        helpers::ZeroOutFOld(latDat);
        helpers::ZeroOutForces(latDat);

        redblood::CellArmy<Kernel> army(*latDat, CellContainer(), cutoff, halo);
        int called = 0;
        auto callback = [cell, &called](std::function<void(CellContainer::value_type)> inserter)
        {
          ++called;
          inserter(cell);
        };
        army.SetCellInsertion(callback);

        army.CallCellInsertion();
        CPPUNIT_ASSERT_EQUAL(1, called);
        CPPUNIT_ASSERT_EQUAL(3ul, army.GetDNC().size());
        CPPUNIT_ASSERT_EQUAL(1ul, army.GetCells().size());
      }

      void CellArmyTests::testFluid2Cell()
      {
        // Checks that positions of cells are updated. Does not check that attendant
        // DNC is.
        auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
        auto const orig = TwoPancakeSamosas<FakeCell>(cutoff);
        auto const normal
          = Facet((*cells.begin())->GetVertices(), (*cells.begin())->GetFacets()[0]).normal();

        LatticePosition gradient;
        Dimensionless non_neg_pop;
        std::function<Dimensionless(PhysicalVelocity const &)> linear, linear_inv;
        std::tie(non_neg_pop, gradient, linear, linear_inv) = helpers::makeLinearProfile(CubeSize(),
                                                                                         latDat,
                                                                                         normal);

        redblood::CellArmy<Kernel> army(*latDat, cells, cutoff, halo);
        army.Fluid2CellInteractions();

        for (size_t i(0); i < cells.size(); ++i)
        {
          auto const disp
            = (*std::next(cells.begin(), i))->GetVertices().front()
            - (*std::next(orig.begin(), i))->GetVertices().front();
          auto i_nodeA = (*std::next(cells.begin(), i))->GetVertices().begin();
          auto i_nodeB = (*std::next(orig.begin(), i))->GetVertices().begin();
          auto const i_end = (*std::next(cells.begin(), i))->GetVertices().end();

          for (; i_nodeA != i_end; ++i_nodeA, ++i_nodeB)
          {
            CPPUNIT_ASSERT(helpers::is_zero( (*i_nodeA - *i_nodeB) - disp));
          }
        }
      }

      void CellArmyTests::testCellOutput()
      {
        auto cell = std::make_shared<FakeCell>(tetrahedron());
        MeshData::Vertices::value_type barycentre;
        typename CellArmy<Kernel>::CellChangeListener callback =
            [&barycentre](const CellContainer & container) {
              barycentre = (*(container.begin()))->GetBarycenter();
        };

        CellContainer intel; intel.insert(cell);
        redblood::CellArmy<Kernel> army(*latDat, intel, cutoff, halo);
        army.AddCellChangeListener(callback);

        army.NotifyCellChangeListeners();
        CPPUNIT_ASSERT_EQUAL(barycentre, cell->GetBarycenter());
      }

      void CellArmyTests::testCellRemoval()
      {
        FlowExtension const outlet(util::Vector3D<Dimensionless>(1, 0, 0),
             LatticePosition(8, 2, 2), 4.0, 4, 1.8);
        auto cell = std::make_shared<FakeCell>(pancakeSamosa());

        helpers::ZeroOutFOld(latDat);
        helpers::ZeroOutForces(latDat);

        CellContainer intel; intel.insert(cell);
        redblood::CellArmy<Kernel> army(*latDat, intel, cutoff, halo);
        army.SetOutlets(std::vector<FlowExtension>(1, outlet));

        // Check status before attempting to remove cell that should *not* be removed
        CPPUNIT_ASSERT_EQUAL(3ul, army.GetDNC().size());
        CPPUNIT_ASSERT_EQUAL(1ul, army.GetCells().size());

        // Check status after attempting to remove cell that should *not* be removed
        army.CellRemoval();
        CPPUNIT_ASSERT_EQUAL(3ul, army.GetDNC().size());
        CPPUNIT_ASSERT_EQUAL(1ul, army.GetCells().size());

        // Check status after attempting to remove cell that should be removed
        *cell += LatticePosition(outlet.origin + outlet.normal * outlet.length * 0.5);
        army.CellRemoval();
        CPPUNIT_ASSERT_EQUAL(0ul, army.GetDNC().size());
        CPPUNIT_ASSERT_EQUAL(0ul, army.GetCells().size());
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellArmyTests);
    }
  }
}

#endif  // ONCE
