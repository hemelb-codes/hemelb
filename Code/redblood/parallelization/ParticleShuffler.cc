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
        std::tuple<bool, proc_t> owner(CellContainer::const_reference cell,
                                       geometry::LatticeData const &latticeData)
        {
          auto const barycenter = cell->GetBarycenter();
          auto const node = latticeData.GetProcIdFromGlobalCoords(barycenter);
          assert(node != BIG_NUMBER2);
          auto const thisrank = latticeData.GetLocalRank();
          return std::make_tuple(bool(thisrank == node), proc_t(node));
        }

        // packs a single cell
        util::Packer& pack(util::Packer& packer, CellContainer::const_reference cell)
        {
          typedef decltype(cell->GetScale()) ScaleType;
          typedef decltype(cell->GetVertices()[0].x) LatticeCoodinateType;
          packer << cell->GetScale();
          packer << cell->GetVertices();
          return packer;
        }

        // unpacks a cell
        util::Packer& unpack(util::Packer& packer, CellContainer::reference cell)
        {
          decltype(cell->GetScale()) scale;
          packer >> scale >> cell->GetVertices();
          cell->SetScale(scale);
          return packer;
        }

        site_t cellSetPackSize(CellContainer const &cells)
        {
          util::Packer::Sizer sizer;
          auto const N = cells.size();
          sizer << N;
          for (auto const& element : cells)
          {
            sizer << element->GetScale() << element->GetVertices();
          }
          return sizer.cast();
        }
        site_t cellSetPackSize(std::map<proc_t, CellContainer> const &cellSet, proc_t node)
        {
          util::Packer::Sizer sizer;
          auto const cells = cellSet.find(node);
          assert(cells != cellSet.end());
          return cellSetPackSize(cells->second);
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
        for (auto const &cell : GetOwnedCells())
        {
          // Check ownership condition
          auto const ownership = GetCellOwner(cell);
          if (std::get<0>(ownership))
            continue;
          // node neighberhood should be declared at start, and not change
          assert(nextCellSwap.count(std::get<1>(ownership)) == 1);
          // if not yet scheduled for moving (anywhere), then add to next shipment
          auto const scheduled =
              std::find_if(currentCellSwap.begin(),
                           currentCellSwap.end(),
                           [&cell](decltype(currentCellSwap)::const_reference aswap)
                           {
                             return aswap.second.count(cell) == 1;
                           });
          if (scheduled == currentCellSwap.end())
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

      util::Packer& ParticleShuffler::Pack(proc_t node, util::Packer& packer) const
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
        packer << cells->second.size();
        for (auto const &cell : cells->second)
        {
          pack(packer, cell);
        }
        return packer;
      }

      util::Packer& ParticleShuffler::Unpack(proc_t node, util::Packer& packer)
      {
#       ifndef NDEBUG
        auto const check = inCommCallOrder.find(node);
        assert(check != inCommCallOrder.end());
        assert(check->second == true);
        check->second = false;
#       endif
        decltype(currentCellSwap.find(0)->second.size()) n;
        packer >> n;
        for (decltype(n) i = 0; i < n; ++i)
        {
          CellContainer::value_type cell = GetEmptyCell();
          unpack(packer, cell);
          AddToOwnedCells(std::move(cell));
        }
        return packer;
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
        for (auto & cells : currentCellSwap)
        {
          for (auto const& cell : cells.second)
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
} // hemelb::redblood::parallel
