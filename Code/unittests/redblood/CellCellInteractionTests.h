//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLCELL_INTERACTION_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLCELL_INTERACTION_TESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"
#include "redblood/CellCell.cc"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellCellInteractionTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (CellCellInteractionTests);
          CPPUNIT_TEST (testBoxHaloTooBig);
          CPPUNIT_TEST (testBoxHalo);
          CPPUNIT_TEST (testAddNodes);
          CPPUNIT_TEST (testAddMeshes);
          CPPUNIT_TEST (testAddCell);
          CPPUNIT_TEST (testRemoveCell);
          CPPUNIT_TEST (testIterator);
          CPPUNIT_TEST (testUpdate);
          CPPUNIT_TEST (testUpdateLentCells);
          CPPUNIT_TEST (testPairIteratorNoPairs);
          CPPUNIT_TEST (testPairIteratorSameMesh);
          CPPUNIT_TEST (testPairIteratorSinglePair);
          CPPUNIT_TEST (testPairIteratorOnePairPerBox);
          CPPUNIT_TEST (testPairIteratorDiagonalBoxes);
          CPPUNIT_TEST (testPairIteratorBoxHalo);CPPUNIT_TEST_SUITE_END();

          LatticeDistance const cutoff = 5.0;
          LatticeDistance const halo = 2.0;

        public:
          void testBoxHaloTooBig();
          void testBoxHalo();
          void testAddNodes();
          void testAddMeshes();
          void testAddCell();
          void testRemoveCell();
          void testIterator();
          void testUpdate();
          void testPairIteratorNoPairs();
          void testPairIteratorSameMesh();
          void testPairIteratorSinglePair();
          void testPairIteratorOnePairPerBox();
          void testPairIteratorDiagonalBoxes();
          void testPairIteratorBoxHalo();
          void testUpdateLentCells();
      };

      void CellCellInteractionTests::testBoxHaloTooBig()
      {
        DivideConquer<int> dnc(cutoff);
        LatticeVector const key(0, 1, -2);
        LatticePosition const position = (LatticePosition(key) + LatticePosition(0.6, 0.5, 0.5))
            * cutoff;

        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(Borders::CENTER)
                                 bitor static_cast<size_t>(Borders::NORTH)
                                 bitor static_cast<size_t>(Borders::SOUTH)
                                 bitor static_cast<size_t>(Borders::EAST)
                                 bitor static_cast<size_t>(Borders::WEST)
                                 bitor static_cast<size_t>(Borders::TOP)
                                 bitor static_cast<size_t>(Borders::BOTTOM),
                             figureNearness(dnc, position, cutoff * 1.001));
        CPPUNIT_ASSERT_EQUAL((size_t) Borders::CENTER bitor (size_t) Borders::TOP,
                             figureNearness(dnc, position, cutoff * 0.499));
      }

      void CellCellInteractionTests::testBoxHalo()
      {
        DivideConquer<int> dnc(5.0);
        LatticeVector const key(0, 1, -2);
        LatticePosition const center = (LatticePosition(key) + LatticePosition(0.5, 0.5, 0.5))
            * dnc.GetBoxSize();
        CPPUNIT_ASSERT(figureNearness(dnc, center, 2.0) == (size_t) Borders::CENTER);

        for (size_t d(1); d < (1 << 6); d <<= 1)
        {
          LatticePosition const disp = direction < LatticePosition::value_type > (d) * 0.6;
          size_t const nearness = figureNearness(dnc, center + disp, 2.0);
          CPPUNIT_ASSERT(nearness == (d bitor (size_t) Borders::CENTER));
        }

        LatticePosition const mult = direction < LatticePosition::value_type
            > (Borders::TOP) * 0.6 + direction < LatticePosition::value_type
            > (Borders::NORTH) * 0.6 + direction < LatticePosition::value_type
            > (Borders::EAST) * 0.6;
        size_t const actual = figureNearness(dnc, center + mult, 2.0);
        size_t const expected = size_t(Borders::TOP) bitor size_t(Borders::NORTH)
            bitor size_t(Borders::EAST) bitor size_t(Borders::CENTER);
        CPPUNIT_ASSERT(actual == expected);
      }

      void CellCellInteractionTests::testAddNodes()
      {
        DivideConquer < CellReference > dnc(cutoff);

        // Adds nodes, last in halo
        LatticePosition const center = LatticePosition(1, 1, 1) * (0.5 * cutoff);
        LatticeDistance const offhalo = (cutoff * 0.5 - 1.1 * halo) / cutoff;
        LatticeDistance const inhalo = (cutoff * 0.5 - 0.9 * halo) / cutoff;
        MeshData::Vertices vertices;
        vertices.push_back(center);
        vertices.push_back(center + LatticePosition(offhalo, 0, 0) * cutoff);
        vertices.push_back(center + LatticePosition(2, offhalo + 3.0, -2) * cutoff);
        vertices.push_back(center + LatticePosition(1, 0, 0) * cutoff + direction
            < LatticePosition::value_type > (Borders::NORTH) * inhalo * cutoff);

        CellContainer cells;
        std::shared_ptr<Cell> intel(new Cell(Mesh(MeshData())));
        cells.insert(intel);

        initializeCells(dnc, vertices, cells.begin(), halo);
        CPPUNIT_ASSERT(dnc.size() == vertices.size());

        DivideConquer<CellReference>::const_range const omega = dnc.equal_range(LatticeVector(0,
                                                                                              0,
                                                                                              0));
        DivideConquer<CellReference>::const_range const empty = dnc.equal_range(LatticeVector(10,
                                                                                              10,
                                                                                              10));
        DivideConquer<CellReference>::const_range const alpha = dnc.equal_range(LatticeVector(2,
                                                                                              3,
                                                                                              -2));
        DivideConquer<CellReference>::const_range const haloed = dnc.equal_range(LatticeVector(1,
                                                                                               0,
                                                                                               0));

        CPPUNIT_ASSERT(std::distance(omega.first, omega.second) == 2);
        CPPUNIT_ASSERT(std::distance(alpha.first, alpha.second) == 1);
        CPPUNIT_ASSERT(std::distance(empty.first, empty.second) == 0);
        CPPUNIT_ASSERT(std::distance(haloed.first, haloed.second) == 1);

        CPPUNIT_ASSERT(helpers::is_zero(omega.first->first));
        CPPUNIT_ASSERT(omega.first->second.cellIterator == cells.begin());
        CPPUNIT_ASSERT(omega.first->second.nodeIndex == 0 or omega.first->second.nodeIndex == 1);
        CPPUNIT_ASSERT(omega.first->second.nearBorder == (size_t) Borders::CENTER);
        DivideConquer<CellReference>::const_iterator other = omega.first;
        ++other;
        CPPUNIT_ASSERT(helpers::is_zero(other->first));
        CPPUNIT_ASSERT(other->second.cellIterator == cells.begin());
        CPPUNIT_ASSERT(other->second.nodeIndex == 0 or other->second.nodeIndex == 1);
        CPPUNIT_ASSERT(other->second.nodeIndex != omega.first->second.nodeIndex);
        CPPUNIT_ASSERT(other->second.nearBorder == (size_t) Borders::CENTER);

        CPPUNIT_ASSERT(alpha.first->first == LatticeVector(2, 3, -2));
        CPPUNIT_ASSERT(alpha.first->second.nodeIndex == 2);
        CPPUNIT_ASSERT(alpha.first->second.cellIterator == cells.begin());
        CPPUNIT_ASSERT(alpha.first->second.nearBorder == (size_t) Borders::CENTER);

        CPPUNIT_ASSERT_EQUAL(LatticeVector::value_type(1), haloed.first->first.x);
        CPPUNIT_ASSERT_EQUAL(LatticeVector::value_type(0), haloed.first->first.y);
        CPPUNIT_ASSERT_EQUAL(LatticeVector::value_type(0), haloed.first->first.z);
        CPPUNIT_ASSERT_EQUAL(site_t(3), haloed.first->second.nodeIndex);
        CPPUNIT_ASSERT(cells.begin() == haloed.first->second.cellIterator);
        CPPUNIT_ASSERT_EQUAL(size_t(Borders::NORTH) bitor size_t(Borders::CENTER),
                             haloed.first->second.nearBorder);
      }

      void checkCell(DivideConquerCells const &dnc, LatticeVector const &key,
                     CellContainer::value_type cell, site_t nbNodes)
      {
        DivideConquerCells::const_range omega = dnc(key);
        CPPUNIT_ASSERT(std::distance(omega.first, omega.second) == nbNodes);

        std::set<int> nodes;

        for (; omega.first != omega.second; ++omega.first)
        {
          CPPUNIT_ASSERT(*omega.first.GetCellReference().cellIterator == cell);
          CPPUNIT_ASSERT(omega.first.GetCellReference().nodeIndex >= 0
              and omega.first.GetCellReference().nodeIndex <= nbNodes);
          CPPUNIT_ASSERT(nodes.count(omega.first.GetCellReference().nodeIndex) == 0);
          nodes.insert(omega.first.GetCellReference().nodeIndex);
        }

        CPPUNIT_ASSERT(nodes.size() == std::size_t(nbNodes));
      }

      void CellCellInteractionTests::testAddMeshes()
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        Mesh pancake = pancakeSamosa();

        CPPUNIT_ASSERT_EQUAL(size_t(2), cells.size());
        auto first = *cells.begin();
        auto second = *std::next(cells.begin());
        if (first->GetBarycenter().GetMagnitude() > second->GetBarycenter().GetMagnitude())
          std::swap(first, second);

        DivideConquerCells dnc(cells, cutoff, halo);
        checkCell(dnc, LatticeVector(0, 0, 0), first, pancake.GetVertices().size());
        checkCell(dnc, LatticeVector(3, 0, 1), second, pancake.GetVertices().size());
      }

      void CellCellInteractionTests::testIterator()
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        DivideConquerCells dnc(cells, cutoff, halo);

        std::set<LatticePosition const *> allnodes;
        typedef MeshData::Vertices::const_iterator vertex_iterator;
        vertex_iterator i_vert = (*cells.begin())->GetVertices().begin();
        vertex_iterator i_vertend = (*cells.begin())->GetVertices().end();

        for (; i_vert != i_vertend; ++i_vert)
        {
          allnodes.insert(& (*i_vert));
        }

        i_vert = (*std::next(cells.begin()))->GetVertices().begin();
        i_vertend = (*std::next(cells.begin()))->GetVertices().end();

        for (; i_vert != i_vertend; ++i_vert)
        {
          allnodes.insert(& (*i_vert));
        }

        DivideConquerCells::const_iterator i_first = dnc.begin();
        DivideConquerCells::const_iterator const i_end = dnc.end();

        for (; i_first != i_end; ++i_first)
        {
          CPPUNIT_ASSERT(not helpers::is_zero(*i_first));
          CPPUNIT_ASSERT(allnodes.count(& (*i_first)) == 1);
          allnodes.erase(& (*i_first));
        }

        CPPUNIT_ASSERT(allnodes.size() == 0);
      }

      void CellCellInteractionTests::testUpdate()
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        DivideConquerCells dnc(cells, cutoff, halo);
        // Figures out which cell is which
        auto firstCell = *cells.begin();
        auto secondCell = *std::next(cells.begin());
        if (firstCell->GetBarycenter().GetMagnitude() > secondCell->GetBarycenter().GetMagnitude())
        {
          std::swap(firstCell, secondCell);
        }

        LatticeVector const zero(0, 0, 0);
        LatticeVector const notzero(3, 0, 1);
        CPPUNIT_ASSERT_EQUAL(dnc.size(), static_cast<unsigned long>(6));
        CPPUNIT_ASSERT_EQUAL(3l, std::distance(dnc(zero).first, dnc(zero).second));
        CPPUNIT_ASSERT_EQUAL(3l, std::distance(dnc(notzero).first, dnc(notzero).second));

        // Move to new box
        LatticeVector const newbox(-1, 1, 2);
        LatticePosition const center(0.5 * cutoff);
        LatticePosition const offhalo = LatticePosition(newbox) * cutoff + center;
        firstCell->GetVertices().front() = offhalo;

        dnc.update();
        CPPUNIT_ASSERT_EQUAL(dnc.size(), static_cast<unsigned long>(6));
        CPPUNIT_ASSERT_EQUAL(3l, std::distance(dnc(notzero).first, dnc(notzero).second));
        CPPUNIT_ASSERT_EQUAL(2l, std::distance(dnc(zero).first, dnc(zero).second));
        CPPUNIT_ASSERT_EQUAL(1l, std::distance(dnc(newbox).first, dnc(newbox).second));
        CPPUNIT_ASSERT(helpers::is_zero(*dnc(newbox).first - offhalo));
        CPPUNIT_ASSERT(not dnc(newbox).first.IsNearBorder());

        // Move near boundary
        LatticePosition const inhalo = LatticePosition(0.25, 0, 0.25) * cutoff + offhalo;
        firstCell->GetVertices().front() = inhalo;
        dnc.update();
        CPPUNIT_ASSERT_EQUAL(dnc.size(), static_cast<unsigned long>(6));
        CPPUNIT_ASSERT_EQUAL(std::distance(dnc(notzero).first, dnc(notzero).second), 3l);
        CPPUNIT_ASSERT_EQUAL(std::distance(dnc(zero).first, dnc(zero).second), 2l);
        CPPUNIT_ASSERT_EQUAL(std::distance(dnc(newbox).first, dnc(newbox).second), 1l);
        CPPUNIT_ASSERT(helpers::is_zero(*dnc(newbox).first - inhalo));
        CPPUNIT_ASSERT(dnc(newbox).first.IsNearBorder());
        CPPUNIT_ASSERT(dnc(newbox).first.IsNearBorder(Borders::TOP));
        CPPUNIT_ASSERT(not dnc(newbox).first.IsNearBorder(Borders::BOTTOM));
        CPPUNIT_ASSERT_EQUAL(dnc(newbox).first.GetNearBorder(),
                             size_t(Borders::TOP) bitor size_t(Borders::EAST)
                                 bitor size_t(Borders::CENTER));
      }

      //! Updating lent cell does not keep extra nodes when number of nodes change in lent cells
      void CellCellInteractionTests::testUpdateLentCells() {
        auto mesh = redblood::pancakeSamosa();
        auto cell = std::make_shared<redblood::Cell>(mesh);
        CellContainer const empty;
        DivideConquerCells dnc(empty, cutoff, halo);

        parallel::ExchangeCells::LentCells lent;
        lent[0].insert(cell);
        dnc.update(parallel::ExchangeCells::ChangedCells(empty, empty, lent));
        CPPUNIT_ASSERT_EQUAL(cell->GetVertices().size(), dnc.size());

        cell->GetVertices().pop_back();
        cell->GetVertices().push_back({4, 4, 4});
        dnc.update(parallel::ExchangeCells::ChangedCells(empty, empty, lent));
        CPPUNIT_ASSERT_EQUAL(cell->GetVertices().size(), dnc.size());

        cell->GetVertices().resize(1);
        dnc.update(parallel::ExchangeCells::ChangedCells(empty, empty, lent));
        CPPUNIT_ASSERT_EQUAL(cell->GetVertices().size(), dnc.size());
      }

      void CellCellInteractionTests::testPairIteratorNoPairs()
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        DivideConquerCells dnc(cells, cutoff, halo);

        // Test when iterating over nothing
        {
          DivideConquerCells::pair_range range(dnc, dnc.end(), dnc.end(), 0.5);
          CPPUNIT_ASSERT(not range.is_valid());
        }

        // Test when no pair are within given range
        {
          DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
          CPPUNIT_ASSERT(not range.is_valid());
        }
      }

      void CellCellInteractionTests::testPairIteratorSameMesh()
      {
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Move one node closer  to the other
        LatticePosition const n0 = (*cells.begin())->GetVertices()[0];
        LatticePosition const n1 = (*cells.begin())->GetVertices()[1];
        (*cells.begin())->GetVertices()[1] = (n1 - n0).GetNormalised() * 0.3 + n0;
        DivideConquerCells dnc(cells, cutoff, halo);

        DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
        DivideConquerCells::const_iterator iterator(dnc.begin());
        CPPUNIT_ASSERT(not range.is_valid());
      }

      void CellCellInteractionTests::testPairIteratorSinglePair()
      {
        // There is only one pair and they are in the same divide and conquer box
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Move one node closer  to the other
        auto const &firstCell = *cells.begin();
        auto const &secondCell = *next(cells.begin());
        LatticePosition const n0 = firstCell->GetVertices()[0];
        LatticePosition const n1 = secondCell->GetVertices()[1];
        secondCell->GetVertices()[1] = (n1 - n0).GetNormalised() * 0.3 + n0;
        DivideConquerCells dnc(cells, cutoff, halo);

        DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
        CPPUNIT_ASSERT(range.is_valid());
        auto const firstNode = range->first.GetCell() == firstCell ?
          *range->first :
          *range->second;
        auto const secondNode = range->first.GetCell() == firstCell ?
          *range->second :
          *range->first;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.x, n0.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.y, n0.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.z, n0.z, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.x, secondCell->GetVertices()[1].x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.y, secondCell->GetVertices()[1].y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.z, secondCell->GetVertices()[1].z, 1e-8);
        CPPUNIT_ASSERT(not ++range);
        CPPUNIT_ASSERT(not range.is_valid());
      }

      void CellCellInteractionTests::testPairIteratorOnePairPerBox()
      {
        // There three pairs and they are each in different boxe, but each contained
        // within one box
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Only one pair, and each in a separate box
        LatticePosition const n0(2 * cutoff - 0.1, 4.5 * cutoff, 4.5 * cutoff);
        LatticePosition const n1(2 * cutoff + 0.1, 4.5 * cutoff, 4.5 * cutoff);
        auto const &firstCell = *cells.begin();
        auto const &secondCell = *next(cells.begin());
        firstCell->GetVertices().front() = n0;
        secondCell->GetVertices().front() = n1;

        DivideConquerCells dnc(cells, cutoff, halo);
        // Checks that fixture is what I think it is
        LatticeVector const N0(1, 4, 4);
        LatticeVector const N1(2, 4, 4);
        CPPUNIT_ASSERT(std::distance(dnc(N0).first, dnc(N0).second) == 1);
        CPPUNIT_ASSERT(std::distance(dnc(N1).first, dnc(N1).second) == 1);

        DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
        CPPUNIT_ASSERT(range.is_valid());

        auto const firstNode = range->first.GetCell() == firstCell ?
          *range->first :
          *range->second;
        auto const secondNode = range->first.GetCell() == firstCell ?
          *range->second :
          *range->first;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.x, n0.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.y, n0.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.z, n0.z, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.x, n1.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.y, n1.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.z, n1.z, 1e-8);
        CPPUNIT_ASSERT(not ++range);
        CPPUNIT_ASSERT(not range.is_valid());
      }

      void CellCellInteractionTests::testPairIteratorDiagonalBoxes()
      {
        // There three pairs and they are each in different boxe, but each contained
        // within one box
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Only one pair, and each in diagonally separate box
        LatticePosition const n0(2 * cutoff - 0.1, 4 * cutoff - 0.1, 8 * cutoff - 0.1);
        LatticePosition const n1(2 * cutoff + 0.1, 4 * cutoff + 0.1, 8 * cutoff + 0.1);
        auto const &firstCell = *cells.begin();
        auto const &secondCell = *next(cells.begin());
        firstCell->GetVertices().front() = n0;
        secondCell->GetVertices().front() = n1;

        DivideConquerCells dnc(cells, cutoff, halo);
        // Checks that fixture is what I think it is
        LatticeVector const N0(1, 3, 7);
        LatticeVector const N1(2, 4, 8);
        CPPUNIT_ASSERT(std::distance(dnc(N0).first, dnc(N0).second) == 1);
        CPPUNIT_ASSERT(std::distance(dnc(N1).first, dnc(N1).second) == 1);
        CPPUNIT_ASSERT_EQUAL(size_t(Borders::CENTER) bitor size_t(Borders::TOP)
                                 bitor size_t(Borders::NORTH) bitor size_t(Borders::EAST),
                             dnc(N0).first.GetNearBorder());
        CPPUNIT_ASSERT_EQUAL(size_t(Borders::CENTER) bitor size_t(Borders::BOTTOM)
                                 bitor size_t(Borders::SOUTH) bitor size_t(Borders::WEST),
                             dnc(N1).first.GetNearBorder());

        DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
        CPPUNIT_ASSERT(range.is_valid());

        auto const firstNode = range->first.GetCell() == firstCell ?
          *range->first :
          *range->second;
        auto const secondNode = range->first.GetCell() == firstCell ?
          *range->second :
          *range->first;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.x, n0.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.y, n0.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.z, n0.z, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.x, n1.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.y, n1.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.z, n1.z, 1e-8);
        CPPUNIT_ASSERT(not ++range);
        CPPUNIT_ASSERT(not range.is_valid());
      }

      void CellCellInteractionTests::testPairIteratorBoxHalo()
      {
        // There three pairs and they are each in different boxe, but each contained
        // within one box
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Only one pair, and each in a separate box
        LatticePosition const n0(2 * cutoff - 0.1, 4.5 * cutoff, 4.5 * cutoff);
        LatticePosition const n1(2 * cutoff + 0.1, 4.5 * cutoff, 4.5 * cutoff);
        auto const &firstCell = *cells.begin();
        auto const &secondCell = *next(cells.begin());
        firstCell->GetVertices().front() = n0;
        secondCell->GetVertices().front() = n1;

        DivideConquerCells dnc(cells, cutoff, halo);
        // Checks that fixture is what I think it is
        LatticeVector const N0(1, 4, 4);
        LatticeVector const N1(2, 4, 4);
        CPPUNIT_ASSERT(std::distance(dnc(N0).first, dnc(N0).second) == 1);
        CPPUNIT_ASSERT(std::distance(dnc(N1).first, dnc(N1).second) == 1);

        DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
        CPPUNIT_ASSERT(range.is_valid());

        auto const firstNode = range->first.GetCell() == firstCell ?
          *range->first :
          *range->second;
        auto const secondNode = range->first.GetCell() == firstCell ?
          *range->second :
          *range->first;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.x, n0.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.y, n0.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(firstNode.z, n0.z, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.x, n1.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.y, n1.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(secondNode.z, n1.z, 1e-8);
        CPPUNIT_ASSERT(not ++range);
        CPPUNIT_ASSERT(not range.is_valid());
      }

      void CellCellInteractionTests::testAddCell()
      {
        auto cell = std::make_shared<Cell>(pancakeSamosa());
        *cell += LatticePosition(1, 1, 1) * cutoff * 0.5;

        CellContainer intelWay;
        auto intelPtr = std::make_shared<Cell>(*cell);
        intelWay.insert(intelPtr);
        DivideConquerCells dnc(intelWay, cutoff, halo);
        const LatticePosition zero(0, 0, 0);
        const LatticePosition one(1, 1, 1);
        CPPUNIT_ASSERT_EQUAL(3ul, dnc.size());
        CPPUNIT_ASSERT_EQUAL(3l, std::distance(dnc(zero).first, dnc(zero).second));

        *cell += LatticePosition(1, 1, 1) * cutoff * 1.5;
        dnc.insert(cell);
        CPPUNIT_ASSERT_EQUAL(6ul, dnc.size());
        CPPUNIT_ASSERT_EQUAL(3l, std::distance(dnc(zero).first, dnc(zero).second));
        CPPUNIT_ASSERT_EQUAL(3l, std::distance(dnc(one).first, dnc(one).second));

        // verifies calling update does nothing
        dnc.update();
        CPPUNIT_ASSERT_EQUAL(6ul, dnc.size());
      }

      void CellCellInteractionTests::testRemoveCell()
      {
        const LatticePosition zero(0, 0, 0);
        auto cells = TwoPancakeSamosas<>(cutoff);
        DivideConquerCells dnc(cells, cutoff, halo);
        auto const zeroBox = dnc(zero);
        CPPUNIT_ASSERT_EQUAL(3l, std::distance(zeroBox.first, zeroBox.second));
        auto const cell = zeroBox.first.GetCell();
        dnc.remove(cell);
        CPPUNIT_ASSERT_EQUAL(3ul, dnc.size());
        CPPUNIT_ASSERT(dnc(zero).first == dnc(zero).second);

        // verifies calling update does nothing
        dnc.update();
        CPPUNIT_ASSERT_EQUAL(3ul, dnc.size());
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellCellInteractionTests);
    }
  }
}

#endif  // ONCE
