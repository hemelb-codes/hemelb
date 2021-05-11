// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_WALLCELLPAIRITERATORTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_WALLCELLPAIRITERATORTESTS_H

#include <iterator>
#include <cppunit/TestFixture.h>
#include "unittests/FourCubeLatticeData.h"
#include "unittests/helpers/HasCommsTestFixture.h"
#include "unittests/helpers/SiteIterator.h"
#include "unittests/helpers/LatticeDataAccess.h"
#include "redblood/WallCellPairIterator.h"
#include "lb/lattices/D3Q15.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class WallCellPairIteratorTests : public helpers::HasCommsTestFixture
      {
          CPPUNIT_TEST_SUITE (WallCellPairIteratorTests);
          CPPUNIT_TEST (testNoWall);
          CPPUNIT_TEST (testOneCellNode);
          CPPUNIT_TEST (testTwoCellNodes);
          CPPUNIT_TEST (testHaloAndNeighboringBoxes);
          CPPUNIT_TEST (testAllWallNodesFound);
          CPPUNIT_TEST (testInteractionPassedOnToFluid<hemelb::redblood::stencil::FourPoint> );CPPUNIT_TEST_SUITE_END();

          LatticeDistance const cutoff = 3.0;
          LatticeDistance const interactionDistance = 0.5;
          LatticeDistance const halo = interactionDistance + 1e-6;
          typedef lb::lattices::D3Q15 Lattice;

        public:
          void setUp()
          {
            latticeData.reset(FourCubeLatticeData::Create(Comms(), 27 + 2));
            for (auto const site : std::cref(*latticeData))
            {
              if (not site.IsWall())
              {
                continue;
              }
              for (Direction d(0); d < Lattice::NUMVECTORS; ++d)
              {
                if (site.HasWall(d))
                {
                  latticeData->SetBoundaryDistance(site.GetIndex(), d, 0.5);
                }
              }
            }
          }

          void testOneCellNode()
          {
            using namespace hemelb::redblood;
            auto const wallDnC = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);

            DivideConquerCells const cellDnC( { cell }, cutoff, interactionDistance);

            auto const last_iterator = std::end(iterate(cellDnC, wallDnC, interactionDistance));
            auto const first_iterator = std::begin(iterate(cellDnC, wallDnC, interactionDistance));
            CPPUNIT_ASSERT(first_iterator != last_iterator);
            CPPUNIT_ASSERT_EQUAL(std::ptrdiff_t(4), std::distance(first_iterator, last_iterator));
            for (auto const item : iterate(cellDnC, wallDnC, interactionDistance))
            {
              CPPUNIT_ASSERT( (item.wallNode - item.cellNode).GetMagnitude() < interactionDistance);
            }
            CPPUNIT_ASSERT_EQUAL(std::ptrdiff_t(0),
                                 std::distance(std::begin(iterate(cellDnC, wallDnC, 0.05)),
                                               std::end(iterate(cellDnC, wallDnC, 0.05))));
          }

          void testNoWall()
          {
            using namespace hemelb::redblood;
            DivideConquer<WallNode> const wallDnC(3e0);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);
            DivideConquerCells const cellDnC( { cell }, cutoff, interactionDistance);

            auto const last_iterator = std::end(iterate(cellDnC, wallDnC, interactionDistance));
            auto const first_iterator = std::begin(iterate(cellDnC, wallDnC, interactionDistance));
            CPPUNIT_ASSERT(first_iterator == last_iterator);
            CPPUNIT_ASSERT_EQUAL(std::ptrdiff_t(0), std::distance(first_iterator, last_iterator));
          }

          void testTwoCellNodes()
          {
            using namespace hemelb::redblood;
            auto const wallDnC = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);
            cell->GetVertices()[1] = LatticePosition(0.49, 3.5 * cutoff, 3.5 * cutoff);

            DivideConquerCells const cellDnC( { cell }, cutoff, interactionDistance);

            CPPUNIT_ASSERT_EQUAL(std::ptrdiff_t(8),
                                 std::distance(std::begin(iterate(cellDnC,
                                                                  wallDnC,
                                                                  interactionDistance)),
                                               std::end(iterate(cellDnC,
                                                                wallDnC,
                                                                interactionDistance))));
          }

          void testHaloAndNeighboringBoxes(Dimensionless where, std::ptrdiff_t howMany)
          {
            using namespace hemelb::redblood;
            auto const wallDnC = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, where * cutoff, where * cutoff);

            DivideConquerCells const cellDnC( { cell }, cutoff, interactionDistance);

            CPPUNIT_ASSERT_EQUAL(howMany,
                                 std::distance(std::begin(iterate(cellDnC,
                                                                  wallDnC,
                                                                  interactionDistance)),
                                               std::end(iterate(cellDnC,
                                                                wallDnC,
                                                                interactionDistance))));
          }

          void testHaloAndNeighboringBoxes()
          {
            // center of a divide and conquer box and in the middle of a square of fluid-sites
            testHaloAndNeighboringBoxes(3.5, 4);
            // center of a divide and conquer box and on top of a fluid-sites
            testHaloAndNeighboringBoxes(3.33333333, 5);
            // // close to range border and in the middle of a square of fluid-sites
            testHaloAndNeighboringBoxes(3.0, 5);
          }

          void testAllWallNodesFound(Dimensionless where)
          {
            using namespace hemelb::redblood;
            auto const wallDnC = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, where * cutoff, where * cutoff);
            DivideConquerCells const cellDnC( { cell }, cutoff, interactionDistance);

            // Creates list of wall-nodes and check unicity
            auto const isSame = [](LatticePosition const &a, LatticePosition const &b)
            {
              return (a - b).GetMagnitude() < 1e-8;
            };
            std::vector<LatticePosition> wallNodes;
            for (auto const nodes : iterate(cellDnC, wallDnC, interactionDistance))
            {
              auto const same = std::find_if(wallNodes.begin(),
                                             wallNodes.end(),
                                             std::bind(isSame,
                                                       nodes.wallNode,
                                                       std::placeholders::_1));
              CPPUNIT_ASSERT(same == wallNodes.end());
              wallNodes.push_back(nodes.wallNode);
            }

            // Loop over all nodes, check whether they are in range, check they are in list
            for (auto const site : std::cref(*latticeData))
            {
              if (not site.IsWall())
              {
                continue;
              }
              for (Direction d(0); d < Lattice::NUMVECTORS; ++d)
              {
                if (site.HasWall(d))
                {
                  auto const node =
                      LatticePosition(Lattice::CX[d], Lattice::CY[d], Lattice::CZ[d]).GetNormalised()
                          * site.GetWallDistance<Lattice>(d) + site.GetGlobalSiteCoords();
                  auto const visited = std::find_if(wallNodes.begin(),
                                                    wallNodes.end(),
                                                    std::bind(isSame, node, std::placeholders::_1));
                  auto const inRange = (node - cell->GetVertices()[0]).GetMagnitude()
                      < interactionDistance;
                  CPPUNIT_ASSERT_EQUAL(visited != wallNodes.end(), inRange);
                  if (visited != wallNodes.end())
                  {
                    wallNodes.erase(visited);
                  }
                }
              }
            }
            // At this point both methods found the same nodes if and only if wallNodes is empty
            CPPUNIT_ASSERT_EQUAL(size_t(0), wallNodes.size());
          }

          void testAllWallNodesFound()
          {
            size_t const N = 10;
            for (size_t i(0); i <= N; ++i)
            {
              testAllWallNodesFound(Dimensionless(i) / Dimensionless(2 * N) + 3.0);
              testAllWallNodesFound(Dimensionless(i) / Dimensionless(2 * N) + 0.0);
            }
          }

          template<class STENCIL> void testInteractionPassedOnToFluid(Dimensionless where)
          {
            using namespace hemelb::redblood;
            auto const wallDnC = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            LatticePosition const node(0.6, where * cutoff, where * cutoff);
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = node;
            DivideConquerCells const cellDnC( { cell }, cutoff, interactionDistance);

            // Set forces to zero
            helpers::ZeroOutForces(static_cast<geometry::LatticeData*>(latticeData.get()));

            // Finds pairs, computes interaction, spread forces to lattice
            addCell2WallInteractions<STENCIL>(DivideConquerCells( { cell }, cutoff, halo),
                                              wallDnC,
                                              Node2NodeForce(1.0, interactionDistance),
                                              *static_cast<geometry::LatticeData*>(latticeData.get()));

            for (auto const site : std::cref(*latticeData))
            {
              auto const d = LatticePosition(site.GetGlobalSiteCoords()) - node;
              CPPUNIT_ASSERT_EQUAL(STENCIL::stencil(d) > 1e-12,
                                   site.GetForce().GetMagnitude() > 1e-12);
            }
          }

          template<class STENCIL> void testInteractionPassedOnToFluid()
          {
            size_t const N = 10;
            for (size_t i(0); i <= N; ++i)
            {
              testInteractionPassedOnToFluid<STENCIL>(Dimensionless(i) / Dimensionless(2 * N)
                  + 3.0);
            }
          }

        private:
          std::unique_ptr<unittests::FourCubeLatticeData> latticeData;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (WallCellPairIteratorTests);
    }
  }
}

#endif  // ONCE
