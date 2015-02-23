//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <iostream>
#include <numeric>
#include "redblood/parallelization/ParticleShuffler.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      namespace
      {
        // id of the proc that should own this cell
        std::tuple<bool, proc_t> owner(
            CellContainer::const_reference cell, geometry::LatticeData const &latticeData)
        {
          auto const barycenter = cell->GetBarycenter();
          auto const node = latticeData.GetProcIdFromGlobalCoords(barycenter);
          assert(node != BIG_NUMBER2);
          auto const thisrank = latticeData.GetLocalRank();
          return {thisrank == node, node};
        }

        // size of the message holding this cell
        site_t cellPackSize(CellContainer::const_reference cell)
        {
          // scale + number of nodes + nodes
          return sizeof(decltype(cell->GetScale()))
            + sizeof(site_t)
            + sizeof(decltype(cell->GetVertices()[0])) * cell->GetNumberOfNodes();
        }

        // packs a single cell
        int8_t* operator<<(int8_t *buffer, CellContainer::const_reference cell)
        {
          typedef decltype(cell->GetScale()) ScaleType;
          typedef decltype(cell->GetVertices()[0].x) LatticeCoodinateType;
          *reinterpret_cast<ScaleType*>(buffer) = cell->GetScale();
          buffer += sizeof(ScaleType);
          *reinterpret_cast<site_t*>(buffer) = static_cast<site_t>(cell->GetNumberOfNodes());
          buffer += sizeof(site_t);
          for(auto const& vertex: cell->GetVertices())
          {
            *reinterpret_cast<LatticeCoodinateType*>(buffer) = vertex.x;
            buffer += sizeof(LatticeCoodinateType);
            *reinterpret_cast<LatticeCoodinateType*>(buffer) = vertex.y;
            buffer += sizeof(LatticeCoodinateType);
            *reinterpret_cast<LatticeCoodinateType*>(buffer) = vertex.z;
            buffer += sizeof(LatticeCoodinateType);
          }
          return buffer;
        }

        // unpacks a cell
        int8_t* operator>>(int8_t *buffer, CellContainer::reference cell)
        {
          typedef decltype(cell->GetScale()) ScaleType;
          typedef decltype(cell->GetVertices()[0].x) LatticeCoodinateType;
          cell->SetScale(*reinterpret_cast<ScaleType*>(buffer));
          buffer += sizeof(ScaleType);
          cell->GetVertices().resize(*reinterpret_cast<site_t*>(buffer));
          buffer += sizeof(site_t);
          for(auto & vertex: cell->GetVertices())
          {
            vertex.x = *reinterpret_cast<LatticeCoodinateType*>(buffer);
            buffer += sizeof(LatticeCoodinateType);
            vertex.y = *reinterpret_cast<LatticeCoodinateType*>(buffer);
            buffer += sizeof(LatticeCoodinateType);
            vertex.z = *reinterpret_cast<LatticeCoodinateType*>(buffer);
            buffer += sizeof(LatticeCoodinateType);
          }
          return buffer;
        }

        site_t cellSetPackSize(std::map<proc_t, CellContainer> const &cellSet, proc_t node)
        {
          auto cSize = [](site_t p, CellContainer::const_reference cell)
          {
            return p + cellPackSize(cell);
          };
          auto const cells = cellSet.find(node);
          assert(cells != cellSet.end());
          return std::accumulate(cells->second.begin(), cells->second.end(), site_t(0), cSize);
        }
      } // anonymous namespace

#     ifndef HEMELB_DOING_UNITTESTS
      ParticleShuffler::ParticleShuffler(std::set<proc_t> const &neighbors)
      {
        for (auto node : neighbors)
        {
          nextCellSwap[node] = CellContainer();
          currentCellSwap[node] = CellContainer();
#         ifndef NDEBUG
            inCommCallOrder[node] = false;
            outCommCallOrder[node] = false;
#         endif
        }
#       ifndef NDEBUG
          callOrder = CallOrder::NONE;
#       endif
      }


      void ParticleShuffler::IdentifyOutOfBounds()
      {
#       ifndef NDEBUG
          assert(callOrder == CallOrder::NONE);
          callOrder = CallOrder::IDENTIFY_CELLS;
#       endif
        for(auto const &cell: GetOwnedCells())
        {
          // Check ownership condition
          auto const ownership = GetCellOwner(cell);
          if(std::get<0>(ownership))
            continue;
          // node neighberhood should be declared at start, and not change
          assert(nextCellSwap.count(std::get<1>(ownership)) == 1);
          // if not yet scheduled for moving (anywhere), then add to next shipment
          auto const scheduled = std::find_if(
              currentCellSwap.begin(), currentCellSwap.end(),
              [&cell](decltype(currentCellSwap)::const_reference aswap)
              {
                return aswap.second.count(cell) == 1;
              }
          );
          if(scheduled == currentCellSwap.end())
            nextCellSwap[std::get<1>(ownership)].insert(cell);
        }
      }

      site_t ParticleShuffler::GetNextCellMessageSize(proc_t node) const
      {
#       ifndef NDEBUG
          assert(callOrder == CallOrder::IDENTIFY_CELLS);
          auto const check = outCommCallOrder.find(node);
          assert(check != outCommCallOrder.end());
          assert(check->second == false);
          check->second = true;
#       endif
        return cellSetPackSize(nextCellSwap, node);
      }

      int8_t* ParticleShuffler::Pack(proc_t node, int8_t* buffer) const
      {
#       ifndef NDEBUG
          assert(callOrder == CallOrder::IDENTIFY_CELLS or callOrder == CallOrder::PACK);
          callOrder = CallOrder::PACK;
          auto const check = outCommCallOrder.find(node);
          assert(check != outCommCallOrder.end());
          assert(check->second == true);
          check->second = false;
#       endif
        auto const cells = currentCellSwap.find(node);
        assert(cells != currentCellSwap.end());
        for(auto const &cell: cells->second)
        {
          buffer = (buffer << cell);
        }
        return buffer;
      }

      int8_t* ParticleShuffler::Unpack(proc_t node, int8_t* buffer)
      {
#       ifndef NDEBUG
          auto const check = inCommCallOrder.find(node);
          assert(check != inCommCallOrder.end());
          assert(check->second == true);
          check->second = false;
#       endif
        auto const bufferSize = incommingCommSize.find(node)->second;
        int8_t * const buffer_end = buffer + bufferSize;
        while(buffer < buffer_end)
        {
          CellContainer::value_type cell = GetEmptyCell();
          buffer = (buffer >> cell);
          AddToOwnedCells(std::move(cell));
        }
        return buffer;
      }

      void ParticleShuffler::SetThisCellMessageSize(proc_t node, site_t messageLength)
      {
#       ifndef NDEBUG
          auto const check = inCommCallOrder.find(node);
          assert(check != inCommCallOrder.end());
          assert(check->second == false);
          check->second = true;
#       endif
          incommingCommSize[node] = messageLength;
      }

      void ParticleShuffler::Next()
      {
#       ifndef NDEBUG
          assert(callOrder == CallOrder::PACK);
          callOrder = CallOrder::NEXT;
#       endif
          // Remove cells that have been sent over
          for(auto & cells: currentCellSwap)
          {
            for(auto const& cell: cells.second)
            {
              RemoveFromOwnedCells(cell);
            }
            cells.second.clear();
          }
          std::swap(nextCellSwap, currentCellSwap);
      }
#     endif
    }
  }
}  // hemelb::redblood::parallel
