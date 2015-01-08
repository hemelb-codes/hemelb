//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLCELL_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLCELL_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"
#include "redblood/CellCell.cc"

namespace hemelb { namespace unittests { namespace redblood {

class CellCellInteractionTests : public CppUnit::TestFixture {

    CPPUNIT_TEST_SUITE(CellCellInteractionTests);
     CPPUNIT_TEST(testBoxHaloTooBig);
     CPPUNIT_TEST(testBoxHalo);
     CPPUNIT_TEST(testAddNodes);
     CPPUNIT_TEST(testAddMeshes);
     CPPUNIT_TEST(testIterator);
     CPPUNIT_TEST(testUpdate);
     CPPUNIT_TEST(testPairIteratorNoPairs);
     CPPUNIT_TEST(testPairIteratorSameMesh);
     CPPUNIT_TEST(testPairIteratorSinglePair);
     CPPUNIT_TEST(testPairIteratorOnePairPerBox);
     CPPUNIT_TEST(testPairIteratorBoxHalo);
    CPPUNIT_TEST_SUITE_END();

    PhysicalDistance const cutoff = 5.0;
    PhysicalDistance const halo = 2.0;

public:
    void testBoxHaloTooBig();
    void testBoxHalo();
    void testAddNodes();
    void testAddMeshes();
    void testIterator();
    void testUpdate();
    void testPairIteratorNoPairs();
    void testPairIteratorSameMesh();
    void testPairIteratorSinglePair();
    void testPairIteratorOnePairPerBox();
    void testPairIteratorBoxHalo();

};

// THIS SPACE IS NECESSARY SINCE UNITTEST CREATES VARIABLES BASED ON LINE #
// THIS SPACE IS NECESSARY SINCE UNITTEST CREATES VARIABLES BASED ON LINE #
// THIS SPACE IS NECESSARY SINCE UNITTEST CREATES VARIABLES BASED ON LINE #
// THIS SPACE IS NECESSARY SINCE UNITTEST CREATES VARIABLES BASED ON LINE #
class CellCellInteractionWithGridTests
  : public helpers::FourCubeBasedTestFixture {
    CPPUNIT_TEST_SUITE(CellCellInteractionWithGridTests);
      CPPUNIT_TEST(testInteraction);
    CPPUNIT_TEST_SUITE_END();

    PhysicalDistance const cutoff = 5.0;
    PhysicalDistance const halo = 2.0;
  public:
    void testInteraction();
  private:
    virtual size_t CubeSize() const { return 32 + 2; }
};

void CellCellInteractionTests :: testBoxHaloTooBig() {
  DivideConquer<int> dnc(cutoff);
  LatticeVector const key(0, 1, -2);
  LatticePosition const position
    = (LatticePosition(key) + LatticePosition(0.6, 0.5, 0.5)) * cutoff;

  CPPUNIT_ASSERT(not figure_nearness(dnc, key, position, cutoff * 0.501));
  CPPUNIT_ASSERT(figure_nearness(dnc, key, position, cutoff * 0.499));
}

void CellCellInteractionTests :: testBoxHalo() {
  DivideConquer<int> dnc(5.0);
  LatticeVector const key(0, 1, -2);
  LatticePosition const center =
    (LatticePosition(key) + LatticePosition(0.5, 0.5, 0.5)) * dnc.GetBoxSize();
  CPPUNIT_ASSERT(figure_nearness(dnc, key, center, 2.0) == 0);
  for(size_t d(1); d < (1 << 6); d <<= 1) {
    LatticePosition const disp = CellReference::directions(d) * 0.6;
    int const nearness = figure_nearness(dnc, key, center + disp, 2.0);
    CPPUNIT_ASSERT(nearness == d);
  }

  LatticePosition const mult =
    CellReference::directions(CellReference::TOP) * 0.6
    + CellReference::directions(CellReference::NORTH) * 0.6
    + CellReference::directions(CellReference::EAST) * 0.6;
  int const actual = figure_nearness(dnc, key, center + mult, 2.0);
  int const expected
    = CellReference::TOP bitor CellReference::NORTH bitor CellReference::EAST;
  CPPUNIT_ASSERT(actual == expected);
}

void CellCellInteractionTests :: testAddNodes() {
  site_t const cellIndex(4);
  DivideConquer<CellReference> dnc(cutoff);

  // Adds nodes, last in halo
  LatticePosition const center = LatticePosition(1, 1, 1) * (0.5 * cutoff);
  PhysicalDistance const offhalo = (cutoff * 0.5 - 1.1 * halo) / cutoff;
  PhysicalDistance const inhalo = (cutoff * 0.5 - 0.9 * halo) / cutoff;
  MeshData::t_Vertices vertices;
  vertices.push_back(center);
  vertices.push_back(center + LatticePosition(offhalo, 0, 0) * cutoff);
  vertices.push_back(center + LatticePosition(2, offhalo + 3.0, -2) * cutoff);
  vertices.push_back(
      center + LatticePosition(1, 0, 0) * cutoff
      + CellReference::directions(CellReference::NORTH) * inhalo * cutoff
  );

  initialize_cells(dnc, vertices, cellIndex, halo);
  CPPUNIT_ASSERT(dnc.size() == vertices.size());

  DivideConquer<CellReference>::const_range const omega
    = dnc.equal_range(LatticeVector(0, 0, 0));
  DivideConquer<CellReference>::const_range const empty
    = dnc.equal_range(LatticeVector(10, 10, 10));
  DivideConquer<CellReference>::const_range const alpha
    = dnc.equal_range(LatticeVector(2, 3, -2));
  DivideConquer<CellReference>::const_range const haloed
    = dnc.equal_range(LatticeVector(1, 0, 0));

  CPPUNIT_ASSERT(std::distance(omega.first, omega.second) == 2);
  CPPUNIT_ASSERT(std::distance(alpha.first, alpha.second) == 1);
  CPPUNIT_ASSERT(std::distance(empty.first, empty.second) == 0);
  CPPUNIT_ASSERT(std::distance(haloed.first, haloed.second) == 1);

  CPPUNIT_ASSERT(helpers::is_zero(omega.first->first));
  CPPUNIT_ASSERT(omega.first->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(
      omega.first->second.nodeIndex == 0 or omega.first->second.nodeIndex == 1
  );
  CPPUNIT_ASSERT(omega.first->second.isNearBorder == 0);
  DivideConquer<CellReference>::const_iterator other = omega.first; ++other;
  CPPUNIT_ASSERT(helpers::is_zero(other->first));
  CPPUNIT_ASSERT(other->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(other->second.nodeIndex == 0 or other->second.nodeIndex == 1);
  CPPUNIT_ASSERT(other->second.nodeIndex != omega.first->second.nodeIndex);
  CPPUNIT_ASSERT(other->second.isNearBorder == 0);

  CPPUNIT_ASSERT(alpha.first->first == LatticeVector(2, 3, -2));
  CPPUNIT_ASSERT(alpha.first->second.nodeIndex == 2);
  CPPUNIT_ASSERT(alpha.first->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(alpha.first->second.isNearBorder == 0);

  CPPUNIT_ASSERT(haloed.first->first == LatticeVector(1, 0, 0));
  CPPUNIT_ASSERT(haloed.first->second.nodeIndex == 3);
  CPPUNIT_ASSERT(haloed.first->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(haloed.first->second.isNearBorder == CellReference::NORTH);
}

void check_cell(
    DivideConquerCells const &_dnc,
    LatticeVector const &_key,
    site_t _cellIndex,
    site_t _nbNodes
) {
  DivideConquerCells::const_range omega = _dnc(_key);
  CPPUNIT_ASSERT(std::distance(omega.first, omega.second) == _nbNodes);

  std::set<int> nodes;
  for(; omega.first != omega.second; ++omega.first) {
    CPPUNIT_ASSERT(omega.first.GetCellReference().cellIndex == _cellIndex);
    CPPUNIT_ASSERT(
        omega.first.GetCellReference().nodeIndex >= 0
        and omega.first.GetCellReference().nodeIndex <= _nbNodes
    );
    CPPUNIT_ASSERT(
        nodes.count(omega.first.GetCellReference().nodeIndex) == 0);
    nodes.insert(omega.first.GetCellReference().nodeIndex);
  }
  CPPUNIT_ASSERT(nodes.size() == _nbNodes);
}

void CellCellInteractionTests :: testAddMeshes() {
  auto cells(TwoPancakeSamosas<>(cutoff));
  Mesh pancake = pancakeSamosa();

  DivideConquerCells dnc(cells, cutoff, halo);
  check_cell(dnc, LatticeVector(0, 0, 0), 0, pancake.GetVertices().size());
  check_cell(dnc, LatticeVector(3, 0, 1), 1, pancake.GetVertices().size());
}

void CellCellInteractionTests :: testIterator() {
  auto cells(TwoPancakeSamosas<>(cutoff));
  DivideConquerCells dnc(cells, cutoff, halo);

  std::set<LatticePosition const*> allnodes;
  typedef MeshData::t_Vertices::const_iterator vertex_iterator;
  vertex_iterator i_vert = cells.front()->GetVertices().begin();
  vertex_iterator i_vertend = cells.front()->GetVertices().end();
  for(; i_vert != i_vertend; ++i_vert)
    allnodes.insert(&(*i_vert));
  i_vert = cells.back()->GetVertices().begin();
  i_vertend = cells.back()->GetVertices().end();
  for(; i_vert != i_vertend; ++i_vert)
    allnodes.insert(&(*i_vert));


  DivideConquerCells::const_iterator i_first = dnc.begin();
  DivideConquerCells::const_iterator const i_end = dnc.end();
  for(; i_first != i_end; ++i_first) {
    CPPUNIT_ASSERT(not helpers::is_zero(*i_first));
    CPPUNIT_ASSERT(allnodes.count(&(*i_first)) == 1);
    allnodes.erase(&(*i_first));
  }
  CPPUNIT_ASSERT(allnodes.size() == 0);
}

void CellCellInteractionTests :: testUpdate() {
  auto cells(TwoPancakeSamosas<>(cutoff));
  DivideConquerCells dnc(cells, cutoff, halo);

  LatticeVector const zero(0, 0, 0);
  LatticeVector const notzero(3, 0, 1);
  CPPUNIT_ASSERT(dnc.size() == 6);
  CPPUNIT_ASSERT(std::distance(dnc(zero).first, dnc(zero).second) == 3);
  CPPUNIT_ASSERT(std::distance(dnc(notzero).first, dnc(notzero).second) == 3);

  // Move to new box
  LatticeVector const newbox(-1, 1, 2);
  LatticePosition const center(0.5 * cutoff);
  LatticePosition const offhalo = LatticePosition(newbox) * cutoff + center;
  cells.front()->GetVertices().front() = offhalo;

  dnc.update();
  CPPUNIT_ASSERT(dnc.size() == 6);
  CPPUNIT_ASSERT(std::distance(dnc(notzero).first, dnc(notzero).second) == 3);
  CPPUNIT_ASSERT(std::distance(dnc(zero).first, dnc(zero).second) == 2);
  CPPUNIT_ASSERT(std::distance(dnc(newbox).first, dnc(newbox).second) == 1);
  CPPUNIT_ASSERT(helpers::is_zero(*dnc(newbox).first - offhalo));
  CPPUNIT_ASSERT(not dnc(newbox).first.IsNearBorder());

  // Move near boundary
  LatticePosition const inhalo
    = LatticePosition(0.25, 0, 0.25) * cutoff + offhalo;
  cells.front()->GetVertices().front() = inhalo;
  dnc.update();
  CPPUNIT_ASSERT(dnc.size() == 6);
  CPPUNIT_ASSERT(std::distance(dnc(notzero).first, dnc(notzero).second) == 3);
  CPPUNIT_ASSERT(std::distance(dnc(zero).first, dnc(zero).second) == 2);
  CPPUNIT_ASSERT(std::distance(dnc(newbox).first, dnc(newbox).second) == 1);
  CPPUNIT_ASSERT(helpers::is_zero(*dnc(newbox).first - inhalo));
  CPPUNIT_ASSERT(dnc(newbox).first.IsNearBorder());
  CPPUNIT_ASSERT(dnc(newbox).first.IsNearBorder(CellReference::TOP));
  CPPUNIT_ASSERT(not dnc(newbox).first.IsNearBorder(CellReference::BOTTOM));
  CPPUNIT_ASSERT(
      dnc(newbox).first.GetNearBorder()
      == (CellReference::TOP bitor CellReference::EAST)
  );
}

void CellCellInteractionTests::testPairIteratorNoPairs() {
  auto cells(TwoPancakeSamosas<>(cutoff));
  DivideConquerCells dnc(cells, cutoff, halo);

  // Test when iterating over nothing
  { DivideConquerCells::pair_range range(dnc, dnc.end(), dnc.end(), 0.5);
    CPPUNIT_ASSERT(not range.is_valid());
  }

  // Test when no pair are within given range
  { DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
    CPPUNIT_ASSERT(not range.is_valid());
  }
}

void CellCellInteractionTests::testPairIteratorSameMesh() {
  auto cells(TwoPancakeSamosas<>(cutoff));

  // Move one node closer  to the other
  LatticePosition const n0 = cells.front()->GetVertices()[0];
  LatticePosition const n1 = cells.front()->GetVertices()[1];
  cells.front()->GetVertices()[1] = (n1 - n0).GetNormalised() * 0.3 + n0;
  DivideConquerCells dnc(cells, cutoff, halo);

  DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
  DivideConquerCells::const_iterator iterator(dnc.begin());
  CPPUNIT_ASSERT(not range.is_valid());
}

void CellCellInteractionTests::testPairIteratorSinglePair() {
  // There is only one pair and they are in the same divide and conquer box
  auto cells(TwoPancakeSamosas<>(cutoff));

  // Move one node closer  to the other
  LatticePosition const n0 = cells.front()->GetVertices()[0];
  LatticePosition const n1 = cells.back()->GetVertices()[1];
  cells.back()->GetVertices()[1] = (n1 - n0).GetNormalised() * 0.3 + n0;
  DivideConquerCells dnc(cells, cutoff, halo);

  DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
  CPPUNIT_ASSERT(range.is_valid());
  CPPUNIT_ASSERT(helpers::is_zero(
        *range->first - cells.front()->GetVertices().front()
  ));
  CPPUNIT_ASSERT(helpers::is_zero(
        *range->second - cells.back()->GetVertices()[1]
  ));
  CPPUNIT_ASSERT(not ++range);
  CPPUNIT_ASSERT(not range.is_valid());
}

void CellCellInteractionTests::testPairIteratorOnePairPerBox() {
  // There three pairs and they are each in different boxe, but each contained
  // within one box
  auto cells(TwoPancakeSamosas<>(cutoff));

  // Only one pair, and each in a separate box
  LatticePosition const n0(2 * cutoff - 0.1, 4.5 * cutoff, 4.5 * cutoff);
  LatticePosition const n1(2 * cutoff + 0.1, 4.5 * cutoff, 4.5 * cutoff);
  cells.front()->GetVertices().front() = n0;
  cells.back()->GetVertices().front() = n1;

  DivideConquerCells dnc(cells, cutoff, halo);
  // Checks that fixture is what I think it is
  LatticeVector const N0(1, 4, 4);
  LatticeVector const N1(2, 4, 4);
  CPPUNIT_ASSERT(std::distance(dnc(N0).first, dnc(N0).second) == 1);
  CPPUNIT_ASSERT(std::distance(dnc(N1).first, dnc(N1).second) == 1);

  DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
  CPPUNIT_ASSERT(range.is_valid());

  CPPUNIT_ASSERT(helpers::is_zero(*range->first - n0));
  CPPUNIT_ASSERT(helpers::is_zero(*range->second - n1));
  CPPUNIT_ASSERT(not ++range);
  CPPUNIT_ASSERT(not range.is_valid());
}

void CellCellInteractionWithGridTests::testInteraction() {
  auto cells = TwoPancakeSamosas<>(cutoff);

  // Place two nodes close enough for interactions
  LatticePosition const n0(15 - 0.1, 15.5, 15.5);
  LatticePosition const n1(15 + 0.1, 15.5, 15.5);
  cells.front()->GetVertices().front() = n0;
  cells.back()->GetVertices().front() = n1;

  // Set forces to zero
  helpers::ZeroOutFOld(latDat);

  // Finds pairs, computes interaction, spread forces to lattice
  addCell2CellInteractions(
      DivideConquerCells(cells, cutoff, halo),
      Node2NodeForce(1.0, halo),
      stencil::FOUR_POINT,
      *latDat
  );

  // By symmetry, there are no forces on the lattice points equidistant from
  // the nodes
  CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 15, 15).GetForce()));
  CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 14, 14).GetForce()));
  CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 16, 16).GetForce()));

  // There are non-zero opposite forces on the following nodes
  CPPUNIT_ASSERT(not helpers::is_zero(latDat->GetSite(14, 15, 15).GetForce()));
  CPPUNIT_ASSERT(helpers::is_zero(
        latDat->GetSite(16, 15, 15).GetForce()
        + latDat->GetSite(14, 15, 15).GetForce()
  ));
  // The forces at (14, 15, 15) should be  in direction (-1, 0, 0)
  CPPUNIT_ASSERT_DOUBLES_EQUAL(
      latDat->GetSite(14, 15, 15).GetForce().Dot(LatticePosition(-1, 0, 0)),
      std::abs(latDat->GetSite(14, 15, 15).GetForce().x),
      1e-8
  );

  // There are non-zero opposite forces on the following nodes
  CPPUNIT_ASSERT(not helpers::is_zero(latDat->GetSite(13, 14, 14).GetForce()));
  CPPUNIT_ASSERT(helpers::is_zero(
        latDat->GetSite(17, 14, 14).GetForce()
        + latDat->GetSite(13, 14, 14).GetForce()
  ));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(
      latDat->GetSite(13, 14, 14).GetForce().Dot(LatticePosition(-1, 0, 0)),
      std::abs(latDat->GetSite(13, 14, 14).GetForce().x),
      1e-8
  );

  // This node is too far away
  CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(12, 15, 15).GetForce()));
}

void CellCellInteractionTests::testPairIteratorBoxHalo() {
  // There three pairs and they are each in different boxe, but each contained
  // within one box
  auto cells = TwoPancakeSamosas<>(cutoff);

  // Only one pair, and each in a separate box
  LatticePosition const n0(2 * cutoff - 0.1, 4.5 * cutoff, 4.5 * cutoff);
  LatticePosition const n1(2 * cutoff + 0.1, 4.5 * cutoff, 4.5 * cutoff);
  cells.front()->GetVertices().front() = n0;
  cells.back()->GetVertices().front() = n1;

  DivideConquerCells dnc(cells, cutoff, halo);
  // Checks that fixture is what I think it is
  LatticeVector const N0(1, 4, 4);
  LatticeVector const N1(2, 4, 4);
  CPPUNIT_ASSERT(std::distance(dnc(N0).first, dnc(N0).second) == 1);
  CPPUNIT_ASSERT(std::distance(dnc(N1).first, dnc(N1).second) == 1);

  DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
  CPPUNIT_ASSERT(range.is_valid());

  CPPUNIT_ASSERT(helpers::is_zero(*range->first - n0));
  CPPUNIT_ASSERT(helpers::is_zero(*range->second - n1));
  CPPUNIT_ASSERT(not ++range);
  CPPUNIT_ASSERT(not range.is_valid());
}


CPPUNIT_TEST_SUITE_REGISTRATION(CellCellInteractionTests);
CPPUNIT_TEST_SUITE_REGISTRATION(CellCellInteractionWithGridTests);
}}}

#endif // ONCE
