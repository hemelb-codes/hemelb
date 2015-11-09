//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <set>
#include <numeric>
#include "util/Iterator.h"
#include "net/MpiError.h"
#include "redblood/parallel/CellParallelization.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      ExchangeCells::ExchangeCells(net::MpiCommunicator const &graphComm) :
        graphComm(graphComm.Duplicate()), step(0)
      {
        auto const neighbors = graphComm.GetNeighbors();
        sendLengths.resize(neighbors.size());
        receiveLengths.resize(neighbors.size());
      }

      void ExchangeCells::PostCellMessageLength(CellParallelization::NodeDistributions const &owned)
      {
        if(step % 3 != 0)
        {
          throw Exception() << "Out-of-order Exchange cell step called\n";
        }
        ++step;
        auto const neighbors = graphComm.GetNeighbors();
        std::fill(sendLengths.begin(), sendLengths.end(), Length{0, 0});
        std::fill(receiveLengths.begin(), receiveLengths.end(), Length{0, 0});
        // Count the number of vertices and cells
        for(auto const & dist: owned)
        {
          for(auto && item: util::zip(neighbors, sendLengths))
          {
            auto const nVertices = dist.second.CountNodes(std::get<0>(item));
            if(nVertices > 0)
            {
              ++std::get<1>(item).nCells;
              std::get<1>(item).nVertices += nVertices;
            }
          }
        }
        // Post message
        HEMELB_MPI_CALL(
          MPI_Ineighbor_alltoall,
          (
            sendLengths.data(), 2, net::MpiDataType<size_t>(),
            receiveLengths.data(), 2, net::MpiDataType<size_t>(),
            graphComm, &lengthRequest
          )
        );
      }

      void ExchangeCells::PostCells() const
      {
      }

      void ExchangeCells::ReceiveCells()
      {
      }

    } // parallel
  } // redblood
}  // hemelb
