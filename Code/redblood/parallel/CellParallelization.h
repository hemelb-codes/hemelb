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
          typedef CellParallelization::NodeDistributions NodeDistributions;
          typedef std::function<int(CellContainer::const_reference)> Ownership;

          ExchangeCells(net::MpiCommunicator const &graphComm)
            : cellCount(graphComm), totalNodeCount(graphComm), nameLengths(graphComm),
              templateNames(graphComm), ownerIDs(graphComm), nodeCount(graphComm),
              cellUUIDs(graphComm), cellScales(graphComm), nodePositions(graphComm)
          {
          }
          //! \brief Computes and posts length of message when sending cells
          //! \param[in] owned: Node distributions of the cells owned by this process
          virtual void PostCellMessageLength(
              NodeDistributions const& distributions, CellContainer const &cells);
          //! \brief Post all owned cells and preps for receiving lent cells
          //! \param[in] distributions tells us for each proc the list of nodes it requires
          //! \param[in] cells a container of cells owned and managed by this process
          //! \param[in] ownership a function to ascertain ownership. It should return the rank of
          //! the owning process in the *world* communicator, according to the position in the cell.
          virtual void PostCells(
              NodeDistributions const &distributions, CellContainer const &cells,
              Ownership const & ownership);
          //! \brief Post all owned cells and preps for receiving lent cells
          //! \param[in] distributions tells us for each proc the list of nodes it requires
          //! \param[in] cells a container of cells owned and managed by this process
          //! \param[in] ownership id of the process owning the cell, corresponding to the rank in
          //! the graph communicator used to build this object.
          virtual void PostCells(
              NodeDistributions const &distributions, CellContainer const &cells,
              std::map<boost::uuids::uuid, proc_t> const & ownership);
          //! Receives messages, reconstructs cells, updates structures and containers
          virtual void ReceiveCells(
              CellContainer &owned, std::map<proc_t, CellContainer> &lent,
              std::shared_ptr<TemplateCellContainer const> const &templateCells);

        protected:
          //! \brief Sends number of cells
          //! \details Using int because it meshes better with sending the number of nodes per cell.
          net::INeighborAllToAll<int> cellCount;
          //! Sends total number of nodes
          net::INeighborAllToAll<int> totalNodeCount;
          //! Send size of message containing template mesh names
          net::INeighborAllToAll<int> nameLengths;
          //! Names of the template meshes
          net::INeighborAllToAllV<char> templateNames;
          //! ID of owner process
          net::INeighborAllToAllV<int> ownerIDs;
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
              NodeDistributions const &distributions, CellContainer const &cells,
              std::map<boost::uuids::uuid, proc_t> const & ownership);
          //! Adds data to local buffers, knowning the neighbor and the cell
          void AddToLocalSendBuffers(
              int neighbor, int nth, int nVertices, proc_t owner,
              CellContainer::const_reference cell,
              NodeCharacterizer::Process2NodesMap::mapped_type const& indices);
          //! Adds local data and send all nodes
          void AddToLocalSendBuffers(
              int neighbor, int nth, proc_t ownerID, CellContainer::const_reference cell);
          //! Adds local data exept nodes
          void AddToLocalSendBuffersAllButNodes(
              int neighbor, int nth, proc_t ownerID, CellContainer::const_reference cell);
          //! Recreates a given cell from messages
          CellContainer::value_type RecreateLentCell(size_t i);
          //! Recreates a given cell from messages
          CellContainer::value_type RecreateOwnedCell(
              size_t i, std::shared_ptr<TemplateCellContainer const> const &templateCells);
      };

    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
