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
          ExchangeCells(net::MpiCommunicator const &graphComm);
          //! \brief Computes and posts length of message when sending cells
          //! \param[in] owned: Node distributions of the cells owned by this process
          virtual void PostCellMessageLength(CellParallelization::NodeDistributions const& owned);
          //! \brief Post all owned cells and preps for receiving lent cells
          virtual void PostCells() const;
          //! Receives cells
          virtual void ReceiveCells();

        protected:
          //! Graph communicator
          net::MpiCommunicator graphComm;
          //! Current step
          int step;
          //! On-going mpi-request
          MPI_Request lengthRequest;

          //! Length message
          struct Length
          {
            //! Number of cells that will be sent
            size_t nCells;
            //! Total number of vertices that will be sent
            size_t nVertices;
          };
          //! Send buffer
          std::vector<Length> sendLengths;
          //! Receive buffer
          std::vector<Length> receiveLengths;
      };
    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
