//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_CELLPARALLELIZATION_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_CELLPARALLELIZATION_H

#include <cppunit/TestFixture.h>
#include <functional>

#include <boost/uuid/uuid_io.hpp>
#include "unittests/redblood/Fixtures.h"
#include "net/MpiCommunicator.h"
#include "redblood/parallel/NodeCharacterizer.h"
#include "redblood/parallel/CellParallelization.h"
#include "util/Iterator.h"
#include <algorithm>

namespace hemelb
{
  namespace unittests
  {
    namespace redblood_parallel
    {
      class ExchangeCells : public redblood::parallel::ExchangeCells
      {
        public:
          ExchangeCells(net::MpiCommunicator const &graphComm) :
              redblood::parallel::ExchangeCells(graphComm)
          {
          }
          ;

#         define HEMELB_MACRO(Name, name, TYPE)   \
            net::TYPE & Get ## Name()             \
            {                                     \
              return name;                        \
            }

          HEMELB_MACRO(CellCount, cellCount, INeighborAllToAll<int>);
          HEMELB_MACRO(TotalNodeCount, totalNodeCount, INeighborAllToAll<int>);
          HEMELB_MACRO(NameLengths, nameLengths, INeighborAllToAll<int>);
          HEMELB_MACRO(NodeCount, nodeCount, INeighborAllToAllV<size_t>);
          HEMELB_MACRO(CellScales, cellScales, INeighborAllToAllV<LatticeDistance>);
          HEMELB_MACRO(OwnerIDs, ownerIDs, INeighborAllToAllV<int>);
          HEMELB_MACRO(CellUUIDs, cellUUIDs, INeighborAllToAllV<unsigned char>);
          HEMELB_MACRO(NodePositions, nodePositions, INeighborAllToAllV<LatticePosition>);
          HEMELB_MACRO(TemplateNames, templateNames, INeighborAllToAllV<char>);
#         undef HEMELB_MACRO
        };

      using namespace hemelb::redblood;
      using namespace hemelb::redblood::parallel;
      class CellParallelizationTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (CellParallelizationTests);
          CPPUNIT_TEST (testCellSwapGetLength);
          CPPUNIT_TEST (testCellSwapPostCells);
          CPPUNIT_TEST (testSingleCellSwapWithRetainedOwnership);
          CPPUNIT_TEST (testSingleCellSwap);
          CPPUNIT_TEST (testUpdateOwnedCells);
          CPPUNIT_TEST (testUpdateNodeDistributions);
          CPPUNIT_TEST (testDistributeCells);
          CPPUNIT_TEST (testLendCells);
          CPPUNIT_TEST (testMotherOfAll);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();
          //! Test message sending number of cells and number of nodes
          void testCellSwapGetLength();
          //! Test messages from swapping cells
          void testCellSwapPostCells();
          //! Test messages from swapping cells while retaining ownership
          void testSingleCellSwapWithRetainedOwnership();
          //! Test messages from swapping cells while retaining ownership
          void testSingleCellSwap();
          //! Checks static function for updating owned cells
          void testUpdateOwnedCells();
          //! Checks static function for updating owned cells
          void testUpdateNodeDistributions();
          //! All cells are in zero at first
          void testDistributeCells();
          //! Two cells are lent from 1 to 0, one of whcih switches ownership
          void testLendCells();
          //! Several cells to and from several processors
          void testMotherOfAll();

          //! Set of nodes affected by given proc
          std::set<proc_t> nodeLocation(LatticePosition const &node);

          //! Center of the region affected by a node
          LatticePosition GetCenter(proc_t rank) const;

          //! Creates map of distributed nodes
          CellParallelization::NodeDistributions GetNodeDistribution(
              CellContainer const &cells) const;

          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(LatticePosition const &pos, boost::uuids::uuid const & uuid,
                                        Dimensionless scale = 1e0, unsigned int depth = 0) const;
          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(LatticePosition const &pos, Dimensionless scale = 1e0,
                                        unsigned int depth = 0) const
          {
            return GetCell(pos, boost::uuids::uuid(), scale, depth);
          }
          std::shared_ptr<Cell> GivenCell(size_t i) const;
          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(Dimensionless scale = 1e0, unsigned int depth = 0) const
          {
            return GetCell(GetCenter(graph ?
              graph.Rank() :
              0),
                           scale, depth);
          }

          //! Id of owning cell
          int Ownership(CellContainer::const_reference cell) const;
          //! Checks that two cells are identical
          void CompareCells(CellContainer::value_type expected,
                            CellContainer::value_type actual) const;
          //! Checks that two cells are identical
          void CompareCells(CellContainer const &expected, CellContainer const & actual) const;
          //! Check that actual contains only vertices from expected affecting this process
          void CompareDistributions(CellContainer::value_type expected,
                                    CellContainer::value_type actual) const;
          void CompareDistributions(CellContainer const & expected,
                                    CellParallelization::LentCells const & actual,
                                    proc_t proc) const;

        protected:
          LatticeDistance const radius = 5;
          net::MpiCommunicator graph;
      };

      void CellParallelizationTests::setUp()
      {
        // setups graph communicator
        auto world = net::MpiCommunicator::World();
        if (world.Size() >= 4)
        {
          std::vector<std::vector<int>> vertices { { 1 }, { 0, 2, 3 }, { 1, 3 }, { 1, 2 } };
          for (int i(4); i < world.Size(); ++i)
          {
            vertices.push_back(std::vector<int> { });
          }
          graph = std::move(world.Graph(vertices));
        }
      }

      //! Get cell centered at given position
      std::shared_ptr<Cell> CellParallelizationTests::GetCell(LatticePosition const &pos,
                                                              boost::uuids::uuid const & uuid,
                                                              Dimensionless scale,
                                                              unsigned int depth) const
      {
        auto cell = std::make_shared<Cell>(icoSphere(depth));
        *cell += pos - cell->GetBarycenter();
        *cell *= scale;
        cell->SetScale(scale);
        cell->SetTag(uuid);
        return cell;
      }

      std::shared_ptr<Cell> CellParallelizationTests::GivenCell(size_t process) const
      {
        size_t const sendto = process > 0 and process < 4 ?
          process - 1 :
          std::numeric_limits<size_t>::max();
        auto const center = GetCenter(sendto);
        auto const scale = 1.0 + 0.1 * static_cast<double>(process);
        boost::uuids::uuid uuid;
        std::fill(uuid.begin(), uuid.end(), static_cast<char>(process));
        auto result = GetCell(center, uuid, scale, process);
        std::ostringstream sstr;
        sstr << process;
        result->SetTemplateName(sstr.str());
        return result;
      }
      ;

      CellParallelization::NodeDistributions CellParallelizationTests::GetNodeDistribution(
          CellContainer const & ownedCells) const
      {
        CellParallelization::NodeDistributions result;
        auto const assessor = std::bind(&CellParallelizationTests::nodeLocation,
                                        *this,
                                        std::placeholders::_1);
        for (auto const&cell : ownedCells)
        {
          result.emplace(std::piecewise_construct,
                         std::forward_as_tuple(cell->GetTag()),
                         std::forward_as_tuple(assessor, cell));
        }
        return result;
      }

      LatticePosition CellParallelizationTests::GetCenter(proc_t rank) const
      {
        switch (rank)
        {
          case 0:
            return
            { 0, 0, 0};
          case 1:
            return
            { 0, 0, 2 * radius};
          case 2:
            return
            { 0, std::sqrt(4 - 1) * radius, 3 * radius};
          case 3:
            return
            { 0, 0, 4 * radius};
          default:
            return
            { 0, 0, static_cast<double>(rank) * radius * 10};
        }
      }

      std::set<proc_t> CellParallelizationTests::nodeLocation(LatticePosition const &node)
      {
        std::set<proc_t> result;
        for (proc_t i(0); i < graph.Size(); ++i)
        {
          if ( (node - GetCenter(i)).GetMagnitude() < radius)
          {
            result.insert(i);
          }
        }
        return result;
      }

      int CellParallelizationTests::Ownership(CellContainer::const_reference cell) const
      {
        std::vector<LatticeDistance> distances;
        auto const ranks = graph ?
          graph.RankMap(net::MpiCommunicator::World()) :
          net::MpiCommunicator::World().RankMap(net::MpiCommunicator::World());
        auto const barycenter = cell->GetBarycenter();
        for (size_t i(0); i < std::size_t(net::MpiCommunicator::World().Size()); ++i)
        {
          auto const d = (GetCenter(ranks.find(i)->second) - barycenter).GetMagnitudeSquared();
          distances.push_back(d);
        };
        return std::min_element(distances.begin(), distances.end()) - distances.begin();
      }

      // Check that cell send each other whole cells
      void CellParallelizationTests::testCellSwapGetLength()
      {
        if (not graph)
        {
          return;
        }
        size_t const sendto = graph.Rank() == 0 ?
                                std::numeric_limits<size_t>::max() :
                              graph.Rank() == 1 ?
                                0 :
                              graph.Rank() == 2 ?
                                1 :
                              graph.Rank() == 3 ?
                                2 :
                                std::numeric_limits<size_t>::max();
        auto const center = GetCenter(sendto);
        CellContainer cells { GetCell(center, 1e0, graph.Rank()) };
        if (graph.Rank() == 0)
        {
          cells.clear();
        }
        auto const dist = GetNodeDistribution(cells);

        ExchangeCells xc(graph);
        CPPUNIT_ASSERT(xc.GetCellCount().GetCommunicator());
        CPPUNIT_ASSERT(xc.GetTotalNodeCount().GetCommunicator());
        auto keepOwnership = [this](CellContainer::const_reference)
        {
          return graph.Rank();
        };
        xc.PostCellMessageLength(dist, cells, keepOwnership);

        // Checks message is correct
        auto const &sendCellCount = xc.GetCellCount().GetSendBuffer();
        auto const &sendTotalNodeCount = xc.GetTotalNodeCount().GetSendBuffer();
        auto const neighbors = graph.GetNeighbors();
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), sendCellCount.size());
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), sendTotalNodeCount.size());
        for (auto const item : util::zip(neighbors, sendCellCount, sendTotalNodeCount))
        {
          auto const sending = std::get<0>(item) == int(sendto);
          size_t const nCells = sending ?
            1 :
            0;
          size_t const nVertices = sending ?
            (*cells.begin())->GetNumberOfNodes() :
            0;
          CPPUNIT_ASSERT_EQUAL(int(nCells), std::get<1>(item));
          CPPUNIT_ASSERT_EQUAL(int(nVertices), int(std::get<2>(item)));
        }

        // Wait for end of request and check received lengths
        xc.GetCellCount().receive();
        xc.GetTotalNodeCount().receive();
        xc.GetNameLengths().receive();
        auto const recvfrom = graph.Rank() == 0 ?
                                1 :
                              graph.Rank() == 1 ?
                                2 :
                              graph.Rank() == 2 ?
                                3 :
                                std::numeric_limits<size_t>::max();

        auto const &receiveCellCount = xc.GetCellCount().GetReceiveBuffer();
        auto const &receiveTotalNodeCount = xc.GetTotalNodeCount().GetReceiveBuffer();
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), receiveCellCount.size());
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), receiveTotalNodeCount.size());
        for (auto const item : util::zip(neighbors, receiveCellCount, receiveTotalNodeCount))
        {
          auto const receiving = std::get<0>(item) == int(recvfrom);
          size_t const nCells = receiving ?
            1 :
            0;
          size_t const nVerts = receiving ?
            GetCell(center, 1e0, recvfrom)->GetNumberOfNodes() :
            0;
          CPPUNIT_ASSERT_EQUAL(int(nCells), std::get<1>(item));
          CPPUNIT_ASSERT_EQUAL(int(nVerts), int(std::get<2>(item)));
        }
      }

      void CellParallelizationTests::testCellSwapPostCells()
      {
        if (not graph)
        {
          return;
        }
        size_t const sendto = graph.Rank() == 0 ?
                                std::numeric_limits<size_t>::max() :
                              graph.Rank() == 1 ?
                                0 :
                              graph.Rank() == 2 ?
                                1 :
                              graph.Rank() == 3 ?
                                2 :
                                std::numeric_limits<size_t>::max();
        auto const center = GetCenter(sendto);
        auto const getScale = [](int process)
        {
          return 1.0 + 0.1 * static_cast<double>(process);
        };
        CellContainer cells { GetCell(center, getScale(graph.Rank()), graph.Rank()) };
        if (graph.Rank() == 0)
        {
          cells.clear();
        }
        auto const dist = GetNodeDistribution(cells);

        ExchangeCells xc(graph);
        auto keepOwnership = [this](CellContainer::const_reference)
        {
          return graph.Rank();
        };
        xc.PostCellMessageLength(dist, cells, keepOwnership);
        xc.PostCells(dist, cells, keepOwnership);

        // check message sizes
        auto const neighbors = graph.GetNeighbors();
        unsigned long const Nsend = graph.Rank() > 0 and graph.Rank() < 4 ?
          1 :
          0;
        unsigned long const Nreceive = graph.Rank() < 3 ?
          1 :
          0;
        unsigned long const uuid_size = sizeof(boost::uuids::uuid);
        CPPUNIT_ASSERT_EQUAL(Nsend, xc.GetCellScales().GetSendBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nreceive, xc.GetCellScales().GetReceiveBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nsend, xc.GetNodeCount().GetSendBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nreceive, xc.GetNodeCount().GetReceiveBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nsend * uuid_size, xc.GetCellUUIDs().GetSendBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nreceive * uuid_size, xc.GetCellUUIDs().GetReceiveBuffer().size());

        if (cells.size() > 0 and graph.Rank() < 4)
        {
          auto const cell = *cells.begin();
          auto const xcScale = xc.GetCellScales().GetSendBuffer()[0];
          auto const xcTag = xc.GetCellUUIDs().GetSendBuffer()[0];
          auto const xcNodes = xc.GetNodeCount().GetSendBuffer()[0];
          CPPUNIT_ASSERT_DOUBLES_EQUAL(cell->GetScale(), xcScale, 1e-8);
          CPPUNIT_ASSERT_EQUAL(*cell->GetTag().begin(), xcTag);
          CPPUNIT_ASSERT_EQUAL(cell->GetNumberOfNodes(), site_t(xcNodes));
        }

        // receive messages
        xc.GetNodeCount().receive();
        xc.GetCellScales().receive();
        xc.GetOwnerIDs().receive();
        xc.GetCellUUIDs().receive();
        xc.GetNodePositions().receive();
        xc.GetTemplateNames().receive();

        if (graph.Rank() < 3)
        {
          auto const scale = getScale(graph.Rank() + 1);
          auto const nNodes = GetCell(center, scale, graph.Rank() + 1)->GetNumberOfNodes();
          CPPUNIT_ASSERT_DOUBLES_EQUAL(scale, xc.GetCellScales().GetReceiveBuffer()[0], 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(nNodes, xc.GetNodeCount().GetReceiveBuffer()[0], 1e-8);
        }
      }

      void CellParallelizationTests::testSingleCellSwapWithRetainedOwnership()
      {
        if (not graph)
        {
          return;
        }

        auto templates =
            std::make_shared<TemplateCellContainer>(TemplateCellContainer { { "1", GivenCell(1) },
                                                                            { "2", GivenCell(2) },
                                                                            { "3", GivenCell(3) } });

        CellContainer owned { GivenCell(graph.Rank()) };
        std::map<proc_t, CellContainer> lent;
        auto const dist = GetNodeDistribution(owned);

        ExchangeCells xc(graph);
        auto keepOwnership = [this](CellContainer::const_reference)
        {
          return graph.Rank();
        };
        xc.PostCellMessageLength(dist, owned, keepOwnership);
        xc.PostCells(dist, owned, keepOwnership);
        auto const result = xc.ReceiveCells(templates);

        if (graph.Rank() < 3)
        {
          auto const & lent = std::get<2>(result);
          CPPUNIT_ASSERT_EQUAL(size_t(1), lent.size());
          CPPUNIT_ASSERT_EQUAL(size_t(1), lent.begin()->second.size());
          CPPUNIT_ASSERT_EQUAL(proc_t(graph.Rank() + 1), lent.begin()->first);
          auto const actual = *lent.begin()->second.begin();
          auto const expected = GivenCell(graph.Rank() + 1);
          CompareCells(actual, expected);
        }
        else
        {
          CPPUNIT_ASSERT_EQUAL(size_t(0), std::get<2>(result).size());
        }
        CPPUNIT_ASSERT_EQUAL(size_t(0), std::get<0>(result).size());
        CPPUNIT_ASSERT_EQUAL(size_t(0), std::get<1>(result).size());
      }

      void CellParallelizationTests::testSingleCellSwap()
      {
        if (not graph)
        {
          return;
        }

        auto templates =
            std::make_shared<TemplateCellContainer>(TemplateCellContainer { { "1", GivenCell(1) },
                                                                            { "2", GivenCell(2) },
                                                                            { "3", GivenCell(3) } });

        CellContainer owned { GivenCell(graph.Rank()) };
        std::map<proc_t, CellContainer> lent;
        auto const dist = GetNodeDistribution(owned);

        ExchangeCells xc(graph);
        auto const ownership = std::bind(&CellParallelizationTests::Ownership,
                                         *this,
                                         std::placeholders::_1);
        xc.PostCellMessageLength(dist, owned, ownership);
        xc.PostCells(dist, owned, ownership);
        auto const result = xc.ReceiveCells(templates);

        CPPUNIT_ASSERT_EQUAL(size_t(0), std::get<2>(result).size());
        if (graph.Rank() < 3)
        {
          CPPUNIT_ASSERT_EQUAL(size_t(1), std::get<0>(result).size());
          auto const actual = *std::get<0>(result).begin();
          auto const expected = GivenCell(graph.Rank() + 1);
          CompareCells(actual, expected);
        }
        if (graph.Rank() < 4 and graph.Rank() > 0)
        {
          CPPUNIT_ASSERT_EQUAL(size_t(1), std::get<1>(result).size());
          CPPUNIT_ASSERT(owned == std::get<1>(result));
        }
      }

      void CellParallelizationTests::testUpdateOwnedCells()
      {
        typedef ExchangeCells::ChangedCells Changes;
        typedef ExchangeCells::LentCells LentCells;
        CellContainer::value_type const cells[4] = { GivenCell(0),
                                                     GivenCell(1),
                                                     GivenCell(2),
                                                     GivenCell(3) };
        CellContainer owned { cells[0], cells[1], cells[2] };

        ExchangeCells::Update(owned,
                              Changes { { cells[3] }, { cells[0], cells[2] }, LentCells { } });
        CPPUNIT_ASSERT( (owned == CellContainer { cells[1], cells[3] }));

        // idem-potent
        ExchangeCells::Update(owned,
                              Changes { { cells[3] }, { cells[0], cells[2] }, LentCells { } });
        CPPUNIT_ASSERT( (owned == CellContainer { cells[1], cells[3] }));
      }

      void CellParallelizationTests::testUpdateNodeDistributions()
      {
        typedef ExchangeCells::ChangedCells Changes;
        typedef ExchangeCells::LentCells LentCells;
        CellContainer::value_type const cells[4] = { GivenCell(0),
                                                     GivenCell(1),
                                                     GivenCell(2),
                                                     GivenCell(3) };
        CellContainer owned { cells[0], cells[1], cells[2] };
        auto distributions = GetNodeDistribution(owned);
        auto const assess = std::bind(&CellParallelizationTests::nodeLocation,
                                      *this,
                                      std::placeholders::_1);

        ExchangeCells::Update(distributions, Changes { { cells[3] },
                                                       { cells[0], cells[2] },
                                                       LentCells { } },
                              assess);
        CPPUNIT_ASSERT_EQUAL(size_t(2), distributions.size());
        CPPUNIT_ASSERT_EQUAL(size_t(1), distributions.count(cells[1]->GetTag()));
        CPPUNIT_ASSERT_EQUAL(size_t(1), distributions.count(cells[3]->GetTag()));

        // idem-potent
        ExchangeCells::Update(distributions, Changes { { cells[3] },
                                                       { cells[0], cells[2] },
                                                       LentCells { } },
                              assess);
        CPPUNIT_ASSERT_EQUAL(size_t(2), distributions.size());
        CPPUNIT_ASSERT_EQUAL(size_t(1), distributions.count(cells[1]->GetTag()));
        CPPUNIT_ASSERT_EQUAL(size_t(1), distributions.count(cells[3]->GetTag()));
      }

      void CellParallelizationTests::testDistributeCells()
      {
        if (not graph)
        {
          return;
        }

        auto const templateCell = GivenCell(0);
        TemplateCellContainer const templates { { templateCell->GetTemplateName(), templateCell } };

        std::vector<int> const sendto = { 0, 0, 1, 2, 3 };
        std::vector<CellContainer::value_type> cells;
        for (size_t i(0); i < sendto.size(); ++i)
        {
          cells.push_back(templateCell->clone());
          *cells[i] += GetCenter(sendto[i]) - cells[i]->GetBarycenter();
          boost::uuids::uuid tag;
          std::fill(tag.begin(), tag.end(), static_cast<unsigned char>(i));
          cells[i]->SetTag(tag);
        }
        *cells[0] += LatticePosition(0.1, 0.1, 0.1);

        auto owned = graph.Rank() == 1 ?
          CellContainer { cells.begin(), cells.end() } :
          CellContainer { };
        auto const dist = GetNodeDistribution(owned);

        ExchangeCells xc(graph);
        auto const ownership = std::bind(&CellParallelizationTests::Ownership,
                                         *this,
                                         std::placeholders::_1);
        xc.PostCellMessageLength(dist, owned, ownership);
        xc.PostCells(dist, owned, ownership);
        auto const result = xc.ReceiveCells(templates);
        xc.Update(owned, result);

        CPPUNIT_ASSERT_EQUAL(size_t(0), std::get<2>(result).size());
        switch (graph.Rank())
        {
          case 0:
            CPPUNIT_ASSERT_EQUAL(size_t(2), owned.size());
            // creating the container should ensure that cells are in the same order
            CompareCells(*CellContainer { cells[0], cells[1] }.begin(), *owned.begin());
            CompareCells(* (++CellContainer { cells[0], cells[1] }.begin()), * (++owned.begin()));
            break;
          case 1:
            CPPUNIT_ASSERT_EQUAL(size_t(1), owned.size());
            CPPUNIT_ASSERT( (owned == CellContainer { cells[2] }));
            break;
          case 2:
            CPPUNIT_ASSERT_EQUAL(size_t(1), owned.size());
            CompareCells(cells[3], *owned.begin());
            break;
          case 3:
            CPPUNIT_ASSERT_EQUAL(size_t(1), owned.size());
            CompareCells(cells[4], *owned.begin());
            break;
          default:
            CPPUNIT_ASSERT_EQUAL(size_t(0), owned.size());
            break;
        }
      }

      void CellParallelizationTests::testLendCells()
      {
        if (not graph)
        {
          return;
        }

        auto const templateCell = GivenCell(0);
        TemplateCellContainer const templates { { templateCell->GetTemplateName(), templateCell } };

        std::vector<int> const sendto = { 0, 1 };
        std::vector<CellContainer::value_type> cells;
        for (size_t i(0); i < sendto.size(); ++i)
        {
          cells.push_back(templateCell->clone());
          // one barycenter on each side of the ownership line
          *cells[i] += (GetCenter(0) + GetCenter(1)) * 0.5 - cells[i]->GetBarycenter();
          *cells[i] += (GetCenter(0) - GetCenter(1)).GetNormalised() * (i != 0 ?
            0.1 :
            -0.1);
          cells[i]->SetScale(1e0 + 0.1 * static_cast<double>(i));
          boost::uuids::uuid tag;
          std::fill(tag.begin(), tag.end(), static_cast<unsigned char>(i));
          cells[i]->SetTag(tag);
        }

        auto owned = graph.Rank() == 1 ?
          CellContainer { cells.begin(), cells.end() } :
          CellContainer { };
        auto const dist = GetNodeDistribution(owned);

        ExchangeCells xc(graph);
        auto const ownership = std::bind(&CellParallelizationTests::Ownership,
                                         *this,
                                         std::placeholders::_1);
        xc.PostCellMessageLength(dist, owned, ownership);
        xc.PostCells(dist, owned, ownership);
        auto const result = xc.ReceiveCells(templates);
        xc.Update(owned, result);

        switch (graph.Rank())
        {
          case 0:
            CPPUNIT_ASSERT_EQUAL(size_t(1), std::get<2>(result).size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), owned.size());
            CompareCells(cells[1], *owned.begin());
            CompareDistributions(cells[0], *std::get<2>(result).begin()->second.begin());
            break;
          case 1:
            CPPUNIT_ASSERT_EQUAL(size_t(1), std::get<2>(result).size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), owned.size());
            CompareCells(cells[0], *owned.begin());
            CompareDistributions(cells[1], *std::get<2>(result).begin()->second.begin());
            break;
          default:
            CPPUNIT_ASSERT_EQUAL(size_t(0), owned.size());
            break;
        }
      }

      void CellParallelizationTests::testMotherOfAll()
      {
        if (not graph)
        {
          return;
        }

        std::vector<CellContainer::value_type> const templateCells = { GivenCell(0),
                                                                       GivenCell(1),
                                                                       GivenCell(2),
                                                                       GivenCell(3) };
        TemplateCellContainer const templates = { { templateCells[0]->GetTemplateName(),
                                                    templateCells[0] },
                                                  { templateCells[1]->GetTemplateName(),
                                                    templateCells[1] },
                                                  { templateCells[2]->GetTemplateName(),
                                                    templateCells[2] },
                                                  { templateCells[3]->GetTemplateName(),
                                                    templateCells[3] } };

        auto const centerOnLine = [=](proc_t i, proc_t j, double alpha)
        {
          auto const a = GetCenter(i), b = GetCenter(j);
          return (b - a) * alpha + a;
        };
        auto const gimmeCell = [=](int n, proc_t i, proc_t j, double alpha, char u)
        {
          auto result = templateCells[n]->clone();
          boost::uuids::uuid tag;
          std::fill(tag.begin(), tag.end(), static_cast<unsigned char>(u));
          result->SetTag(tag);
          *result += centerOnLine(i, j, alpha) - result->GetBarycenter();
          return CellContainer::value_type(std::move(result));
        };
        std::vector<CellContainer::value_type> const cells {
        // owned by 0
                                                             gimmeCell(0, 0, 1, 0.45, 0),
                                                             gimmeCell(0, 0, 1, 0.48, 1),
                                                             gimmeCell(1, 0, 1, 0.55, 2),
                                                             // owned by 1
                                                             gimmeCell(1, 1, 0, 0.0, 3),
                                                             gimmeCell(1, 1, 2, 0.48, 4),
                                                             gimmeCell(2, 1, 2, 0.55, 5),
                                                             gimmeCell(3, 1, 3, 0.48, 6),
                                                             gimmeCell(2, 1, 3, 0.54, 7),
                                                             gimmeCell(3, 1, 3, 0.49, 8),
                                                             gimmeCell(0, 1, 3, 0.54, 9),
                                                             gimmeCell(3, 1, 3, 0.0, 10),
                                                             gimmeCell(0, 1, 3, 1.0, 11),
                                                             // owned by 2
                                                             gimmeCell(2, 2, 1, 0.48, 12),
                                                             gimmeCell(1, 2, 1, 0.47, 13),
                                                             gimmeCell(2, 2, 3, 0.45, 14) };

        auto owned = graph.Rank() == 0 ?
                       CellContainer { cells.begin(), cells.begin() + 3 } :
                     graph.Rank() == 1 ?
                       CellContainer { cells.begin() + 3, cells.begin() + 3 + 9 } :
                     graph.Rank() == 2 ?
                       CellContainer { cells.begin() + 3 + 9, cells.end() } :
                       CellContainer { };
        auto const dist = GetNodeDistribution(owned);

        ExchangeCells xc(graph);
        auto const ownership = std::bind(&CellParallelizationTests::Ownership,
                                         *this,
                                         std::placeholders::_1);
        xc.PostCellMessageLength(dist, owned, ownership);
        xc.PostCells(dist, owned, ownership);
        auto const result = xc.ReceiveCells(templates);
        xc.Update(owned, result);

        switch (graph.Rank())
        {
          case 0:
            CompareCells( { cells[0], cells[1] }, owned);
            CompareDistributions( { cells[2] }, std::get<2>(result), 1);
            break;
          case 1:
            CompareCells( { cells[2], cells[3], cells[4], cells[6], cells[8], cells[10] }, owned);
            CompareDistributions( { cells[0], cells[1] }, std::get<2>(result), 0);
            CompareDistributions( { cells[5], cells[12], cells[13] }, std::get<2>(result), 2);
            CompareDistributions( { cells[9] }, std::get<2>(result), 3);
            break;
          case 2:
            CompareCells( { cells[5], cells[12], cells[13], cells[14] }, owned);
            CompareDistributions( { cells[4] }, std::get<2>(result), 1);
            CPPUNIT_ASSERT_EQUAL(size_t(0), std::get<2>(result).count(3));
            break;
          case 3:
            CompareCells( { cells[7], cells[9], cells[11] }, owned);
            CompareDistributions( { cells[6], cells[8] }, std::get<2>(result), 1);
            CompareDistributions( { cells[14] }, std::get<2>(result), 2);
            break;
          default:
            CPPUNIT_ASSERT_EQUAL(size_t(0), owned.size());
            break;
        }
      }

      void CellParallelizationTests::CompareCells(CellContainer::value_type expected,
                                                  CellContainer::value_type actual) const
      {
        CPPUNIT_ASSERT_EQUAL(expected->GetTag(), actual->GetTag());
        CPPUNIT_ASSERT_EQUAL(expected->GetTemplateName(), actual->GetTemplateName());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected->GetScale(), actual->GetScale(), 1e-8);
        CPPUNIT_ASSERT_EQUAL(expected->GetNumberOfNodes(), actual->GetNumberOfNodes());
        for (auto item : util::zip(expected->GetVertices(), actual->GetVertices()))
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
        }
      }
      void CellParallelizationTests::CompareCells(CellContainer const &expected,
                                                  CellContainer const & actual) const
      {
        CPPUNIT_ASSERT_EQUAL(expected.size(), actual.size());
        for (auto const &item : util::zip(expected, actual))
        {
          CompareCells(std::get<0>(item), std::get<1>(item));
        }
      }

      void CellParallelizationTests::CompareDistributions(
          CellContainer::value_type expectedCell, CellContainer::value_type actualCell) const
      {
        CPPUNIT_ASSERT_EQUAL(expectedCell->GetTag(), actualCell->GetTag());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedCell->GetScale(), actualCell->GetScale(), 1e-8);
        CPPUNIT_ASSERT_EQUAL(expectedCell->GetTemplateName(), actualCell->GetTemplateName());
        auto const distributions = GetNodeDistribution( { expectedCell });
        auto const distribution = distributions.find(expectedCell->GetTag())->second;
        auto const expectedNodes = distribution.CountNodes(graph.Rank()) > 0 ?
          distribution[graph.Rank()] :
          decltype(distribution[graph.Rank()]) { };
        CPPUNIT_ASSERT_EQUAL(size_t(expectedNodes.size()), size_t(actualCell->GetNumberOfNodes()));
        auto const &expectedV = expectedCell->GetVertices();
        auto const &actualV = actualCell->GetVertices();
        for (auto const item : util::enumerate(expectedNodes))
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedV[item.value].x, actualV[item.index].x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedV[item.value].y, actualV[item.index].y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedV[item.value].z, actualV[item.index].z, 1e-8);
        };
      }

      void CellParallelizationTests::CompareDistributions(
          CellContainer const & expected, CellParallelization::LentCells const & actual,
          proc_t proc) const
      {
        auto const i_proc = actual.find(proc);
        CPPUNIT_ASSERT(i_proc != actual.end());
        CPPUNIT_ASSERT_EQUAL(expected.size(), i_proc->second.size());
        for (auto const &item : util::zip(expected, i_proc->second))
        {
          CompareDistributions(std::get<0>(item), std::get<1>(item));
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellParallelizationTests);
    }
  }
}

#endif  // ONCE
