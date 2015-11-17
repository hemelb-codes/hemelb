//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_REDBLOODE_PARALLELIZATION_CELL_PARALLELIZATION
#define HEMELB_REDBLOODE_PARALLELIZATION_CELL_PARALLELIZATION

#include <boost/uuid/uuid.hpp>
#include <map>

#include "redblood/parallel/NodeCharacterizer.h"
#include "redblood/Cell.h"
#include "net/MpiCommunicator.h"
#include "net/INeighborAllToAll.h"
#include "net/INeighborAllToAllV.h"


namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      class CellParallelization
      {
        public:
          //! Type of the object holding distributions
          typedef std::map<boost::uuids::uuid, NodeCharacterizer> NodeDistributions;
          //! Creates a cell-parallelization object
          CellParallelization(net::MpiCommunicator const &comm)
            : comm(comm.Duplicate())
          {
          }
          //! Creates a cell-parallelization object
          CellParallelization(net::MpiCommunicator &&comm)
            : comm(std::move(comm))
          {
          }
          //! Adds an owned cell
          void AddCell(CellContainer::const_reference cell)
          {
            owned.insert(cell);
          }

        protected:
          //! Graph communicator defining neighberhood over which cells can be owned
          net::MpiCommunicator comm;
          //! Cells owned by this process
          CellContainer owned;
          //! Distribution of the cells owned by this process
          NodeDistributions distributions;
          //! Cells lent by other proces
          std::map<proc_t, CellContainer> lent;
      };

      class ExchangeCells
      {
        public:
          ExchangeCells(net::MpiCommunicator const &graphComm)
            : cellCount(graphComm), totalNodeCount(graphComm), nodeCount(graphComm),
              cellUUIDs(graphComm), cellScales(graphComm), nodePositions(graphComm)
          {
          }
          //! \brief Computes and posts length of message when sending cells
          //! \param[in] owned: Node distributions of the cells owned by this process
          virtual void PostCellMessageLength(CellParallelization::NodeDistributions const& owned);
          //! \brief Post all owned cells and preps for receiving lent cells
          virtual void PostCells(
            CellParallelization::NodeDistributions const &owned, CellContainer const &cells);
          //! Receives cells
          virtual void ReceiveCells();

        protected:
          //! \brief Sends number of cells
          //! \details Using int because it meshes better with sending the number of nodes per cell.
          net::INeighborAllToAll<int> cellCount;
          //! Sends total number of nodes
          net::INeighborAllToAll<size_t> totalNodeCount;
          //! Sends number of nodes per cell
          net::INeighborAllToAllV<size_t> nodeCount;
          //! Cell uuids
          net::INeighborAllToAllV<unsigned char> cellUUIDs;
          //! Sends scales
          net::INeighborAllToAllV<LatticeDistance> cellScales;
          //! Sends positions
          net::INeighborAllToAllV<LatticeDistance> nodePositions;

          //! Number of nodes to send to each neighboring process
          void SetupLocalSendBuffers(
              CellParallelization::NodeDistributions const &owned, CellContainer const &cells);
          void AddToLocalSendBuffers(
              int neighbor, int nth, int nVertices, CellContainer::const_reference cell);
      };

    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
