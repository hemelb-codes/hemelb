// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.


#include <functional>
#include <catch2/catch.hpp>

#include "net/MpiCommunicator.h"
#include "net/INeighborAllToAll.h"
#include "net/INeighborAllToAllV.h"
#include "util/Iterator.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::net;

    TEST_CASE("NeighborCommTests") {
      net::MpiCommunicator graph;

      // setups graph communicator
      auto world = net::MpiCommunicator::World();
      if(world.Size() >= 4) {
	std::vector<std::vector<int>> vertices{{1}, {0, 2, 3}, {1, 3}, {1, 2}};
	for(int i = 4; i < world.Size(); ++i) {
	  vertices.push_back(std::vector<int>{});
	}
	graph = world.Graph(vertices);
      }

      SECTION("Sends one integer per proc") {
        if(not graph)
        {
          return;
        }

        // Just a message such what i sends to j is different from what j sends to i
        auto const message = [](int proci, int procj)
        {
          return proci + procj + (proci > procj ? 10: proci < procj ? -10: 0);
        };

        // Prepares buffer. One message per object.
        auto const neighbors = graph.GetNeighbors();
        INeighborAllToAll<int> all2all(graph);
        all2all.GetSendBuffer().resize(neighbors.size());
        for(auto const item: util::zip(all2all.GetSendBuffer(), neighbors))
        {
          std::get<0>(item) = message(graph.Rank(), std::get<1>(item));
        }
        // Send message
        all2all.send();
        // Wait for end of request and check received lengths
        all2all.receive();

	REQUIRE(neighbors.size() == all2all.GetReceiveBuffer().size());
        for(auto const item: util::zip(all2all.GetReceiveBuffer(), neighbors))
        {
          auto const received = std::get<0>(item);
          auto const neighbor = std::get<1>(item);
	  REQUIRE(message(neighbor, graph.Rank()) == received);
        }
      }

      SECTION("Sends several doubles per proc") {
        if(not graph)
        {
          return;
        }
        // Just a message such what i sends to j is different from what j sends to i
        auto const message = [](int proci, int procj, int n = 0)
        {
          return n + proci + procj + (proci > procj ? 100: proci < procj ? -100: 0);
        };

        // Prepares buffer. One message per object.
        auto const neighbors = graph.GetNeighbors();
        INeighborAllToAll<int> all2all(graph);
        for(auto const n: neighbors)
        {
          auto const r = graph.Rank();
          all2all.AddToBuffer(n, {message(r, n), message(r, n, 1), message(r, n, 2)});
        }
        // Send message
        all2all.send();
        // Wait for end of request and check received lengths
        all2all.receive();

        auto const & received = all2all.GetReceiveBuffer();
        REQUIRE(neighbors.size() * 3 == received.size());
        for(auto const item: util::enumerate(neighbors))
        {
          REQUIRE(message(item.value, graph.Rank()) == received[item.index * 3]);
          REQUIRE(message(item.value, graph.Rank(), 1) == received[item.index * 3 + 1]);
          REQUIRE(message(item.value, graph.Rank(), 2) == received[item.index * 3 + 2]);
        }
      }

      SECTION("testVariableInts") {
	// Sends different number of ints from each proc
	// 
	// Each proc expects to know how many it is gettin from each before hand
        if(not graph)
        {
          return;
        }
        // Just a message such what i sends to j is different from what j sends to i
        auto const message = [](int proci, int procj, int n = 0)
        {
          return n + proci + procj + (proci > procj ? 100: proci < procj ? -100: 0);
        };
        // If X sends to Y, then how many it does send
        auto const XSendsToY = [](int proci, int procj)
        {
          return proci >= 4 or procj >= 4 ? 0:
              proci > procj ? proci + procj: proci + procj + 10;
        };

        // Prepares arrays how many elements are sent from each to each
        auto const neighbors = graph.GetNeighbors();
        std::vector<int> sendCounts, recvCounts;
        int nReceived = 0;
        for(auto neighbor: neighbors)
        {
          sendCounts.push_back(XSendsToY(graph.Rank(), neighbor));
          recvCounts.push_back(XSendsToY(neighbor, graph.Rank()));
          nReceived += recvCounts.back();
        }

        // Prepares array with elements to send
        INeighborAllToAllV<int> all2all(graph, sendCounts, recvCounts);
        for(auto const neighbor: neighbors)
        {
          std::vector<int> stuff;
          for(int i(0); i < XSendsToY(graph.Rank(), neighbor); ++i)
          {
            stuff.push_back(message(graph.Rank(), neighbor, i));
          }
          all2all.insertSend(neighbor, stuff);
        }

        // Send message
        all2all.send();
        // Wait for end of request and check received lengths
        all2all.receive();

        // Check result
        REQUIRE(size_t(nReceived) == all2all.GetReceiveBuffer().size());
        auto received = all2all.GetReceiveBuffer().begin();
        for(auto const item: util::zip(neighbors, recvCounts))
        {
          for(int i(0); i < std::get<1>(item); ++i, ++received)
          {
	    REQUIRE(message(std::get<0>(item), graph.Rank(), i) == *received);
          }
        }

      }

      SECTION("Tests mapping of ranks from one comm to another") {

        if(not graph)
        {
          return;
        }

        int const N = graph.Size();
        auto newRank = [N](int i) { return (i + 5) % N; };
        auto const newComm = graph.Split(0, newRank(graph.Rank()));
        auto const rankMap = newComm.RankMap(graph);
        for(int i(0); i < graph.Size(); ++i)
        {
          REQUIRE(size_t(1) == rankMap.count(newRank(i)));
          REQUIRE(i == rankMap.find(newRank(i))->second);
        }
      }      
    }

  }
}
