//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_PARALLEL_PARTICLESHUFFLER_H
#define HEMELB_REDBLOOD_PARALLEL_PARTICLESHUFFLER_H

#include <tuple>

#include "geometry/LatticeData.h"
#include "redblood/Cell.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      //! Shuffles particles from one to process to another
      class ParticleShuffler
      {
        public:
          //! Object that determines how to shuffle particles
          //! \param[in] neighbors: Defines nodes this proc will be in communication with
          ParticleShuffler(std::set<proc_t> const &neighbors);
          virtual ~ParticleShuffler()
          {
          }

          //! Sets cells in nextCellSwap
          void IdentifyOutOfBounds();
          //! Defines (lengths of) next cell-message
          site_t GetNextCellMessageSize(proc_t node) const;
          //! Set size of next incoming message
          void SetThisCellMessageSize(proc_t node, site_t messageLength);
          //! Pack's current cell message
          //! \param[in] node: node for which to pack a message
          //! \param[in] buffer: must have correct size, given by call to GetNextCellMessageSize for
          //! given proc.
          int8_t* Pack(proc_t node, int8_t *buffer) const;
          //! Unpack's cell message received from given proc
          //!
          //! It is assumed that the size of the message to unpack has been set previously by a call
          //! to SetCellMessageSize(currentCellSwap, node);
          //!
          //! \param[in] node: node from which to unpack a message
          //! \param[in] buffer: must have correct size, given by call to GetNextCellMessageSize for
          //!            given proc.
          int8_t* Unpack(proc_t node, int8_t* buffer);
          //! Prepares for next iteration
          void Next();


        protected:
          //! Figures out who owns this object
          //! \returns tuple (self owned, rank). The first item indicates whether the cell is still
          //! owned by this object. If not the second item indicates the rank of the owner.
          virtual std::tuple<bool, proc_t> GetCellOwner(CellContainer::const_reference) const = 0;
          //! Const reference to cells owned by this proc
          virtual CellContainer const& GetOwnedCells() const = 0;
          //! Returns an empty cell that will be filled by unpacking a message from another proc
          virtual CellContainer::value_type GetEmptyCell() const = 0;
          //! Adds a cell as owned
          virtual void AddToOwnedCells(CellContainer::value_type&&) = 0;
          //! Removes owned cell
          virtual void RemoveFromOwnedCells(CellContainer::const_reference) = 0;

        protected:
          //! Holds cells to swap in next LB iteration
          std::map<proc_t, CellContainer> nextCellSwap;
          //! Holds cells to swap in current iteration
          std::map<proc_t, CellContainer> currentCellSwap;
          //! Holds size of incomming iteration
          std::map<proc_t, size_t> incommingCommSize;

        private:
#         ifndef NDEBUG
            //! True if incomming size has been set
            std::map<proc_t, bool> inCommCallOrder;
            //! True if outgoing size has been checked
            mutable std::map<proc_t, bool> outCommCallOrder;
            //! True if IdentifyOutOfBounds called in right order
            enum class CallOrder
            {
              NONE,
              IDENTIFY_CELLS,
              PACK,
              NEXT = NONE
            } mutable callOrder;
#         endif
      };
    }
  }
}
#endif
