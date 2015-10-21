//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_WALL_NODE_DNC_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_WALL_NODE_DNC_TESTS_H

#include <cppunit/TestFixture.h>
#include "unittests/FourCubeLatticeData.h"
#include "unittests/helpers/HasCommsTestFixture.h"
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
          CPPUNIT_TEST (testOneCellNode);
          CPPUNIT_TEST (testTwoCellNodes);
          CPPUNIT_TEST (testHaloAndNeighboringBoxes);
          CPPUNIT_TEST (testAllWallNodesFound);
          CPPUNIT_TEST_SUITE_END();

          LatticeDistance const cutoff = 3.0;
          LatticeDistance const interactionDistance = 0.5;
          LatticeDistance const halo = interactionDistance + 1e-6;
          typedef lb::lattices::D3Q15 Lattice;

        public:
          void setUp()
          {
            latticeData.reset(FourCubeLatticeData::Create(Comms(), 27+2));
            for(site_t i(0); i < latticeData->GetLocalFluidSiteCount(); ++i)
            {
              auto const site = latticeData->GetSite(i);
              if(not site.IsWall())
              {
                continue;
              }
              for(Direction d(0); d < Lattice::NUMVECTORS; ++d)
              {
                if(site.HasWall(d))
                {
                  latticeData->SetBoundaryDistance(i, d, 0.5);
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
            *cell += 100;
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);

            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            auto const last_iterator = std::end(iterate(cellDnC, wallDnC, interactionDistance));
            auto const first_iterator = std::begin(iterate(cellDnC, wallDnC, interactionDistance));
            CPPUNIT_ASSERT(first_iterator != last_iterator);
            CPPUNIT_ASSERT_EQUAL(std::ptrdiff_t(4), std::distance(first_iterator, last_iterator));
            for(auto const item: iterate(cellDnC, wallDnC, interactionDistance))
            {
              CPPUNIT_ASSERT((item.wallNode - item.cellNode).GetMagnitude() < interactionDistance);
            }
            CPPUNIT_ASSERT_EQUAL(
                std::ptrdiff_t(0),
                std::distance(
                  std::begin(iterate(cellDnC, wallDnC, 0.05)),
                  std::end(iterate(cellDnC, wallDnC, 0.05))
                )
            );
          }

          void testTwoCellNodes()
          {
            using namespace hemelb::redblood;
            auto const wallDnC = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += 100;
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);
            cell->GetVertices()[1] = LatticePosition(0.49, 3.5 * cutoff, 3.5 * cutoff);

            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            CPPUNIT_ASSERT_EQUAL(
                std::ptrdiff_t(8),
                std::distance(
                  std::begin(iterate(cellDnC, wallDnC, interactionDistance)),
                  std::end(iterate(cellDnC, wallDnC, interactionDistance))
                )
            );
          }

          void testHaloAndNeighboringBoxes(Dimensionless where, std::ptrdiff_t howMany)
          {
            using namespace hemelb::redblood;
            auto const wallDnC = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += 100;
            cell->GetVertices()[0] = LatticePosition(0.5, where * cutoff, where * cutoff);

            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            CPPUNIT_ASSERT_EQUAL(
                howMany,
                std::distance(
                  std::begin(iterate(cellDnC, wallDnC, interactionDistance)),
                  std::end(iterate(cellDnC, wallDnC, interactionDistance))
                )
            );
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
            *cell += 100;
            cell->GetVertices()[0] = LatticePosition(0.5, where * cutoff, where * cutoff);
            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            // Creates list of wall-nodes and check unicity
            auto const isSame = [](LatticePosition const &a, LatticePosition const &b)
            {
              return (a - b).GetMagnitude() < 1e-8;
            };
            std::vector<LatticePosition> wallNodes;
            for(auto const nodes: iterate(cellDnC, wallDnC, interactionDistance))
            {
              auto const same = std::find_if(
                  wallNodes.begin(), wallNodes.end(),
                  std::bind(isSame, nodes.wallNode, std::placeholders::_1)
              );
              CPPUNIT_ASSERT(same == wallNodes.end());
              wallNodes.push_back(nodes.wallNode);
            }

            // Loop over all nodes, check whether they are in range, check they are in list
            for(site_t i(0); i < latticeData->GetLocalFluidSiteCount(); ++i)
            {
              auto const site = latticeData->GetSite(i);
              if(not site.IsWall())
              {
                continue;
              }
              for(Direction d(0); d < Lattice::NUMVECTORS; ++d)
              {
                if(site.HasWall(d))
                {
                  auto  const node
                    = LatticePosition(Lattice::CX[d], Lattice::CY[d], Lattice::CZ[d])
                        .GetNormalised() * site.GetWallDistance<Lattice>(d)
                      + site.GetGlobalSiteCoords();
                  auto const visited = std::find_if(
                      wallNodes.begin(), wallNodes.end(),
                      std::bind(isSame, node, std::placeholders::_1)
                  );
                  auto const inRange
                    = (node - cell->GetVertices()[0]).GetMagnitude() < interactionDistance;
                  CPPUNIT_ASSERT_EQUAL(visited != wallNodes.end(), inRange);
                  if(visited != wallNodes.end())
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
            for(size_t i(0); i <= N; ++i)
            {
              testAllWallNodesFound(Dimensionless(i)/Dimensionless(2*N) + 3.0);
              testAllWallNodesFound(Dimensionless(i)/Dimensionless(2*N) + 0.0);
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
