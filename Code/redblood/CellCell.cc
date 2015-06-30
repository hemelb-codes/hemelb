//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <vector>
#include "redblood/CellCell.h"
#include "redblood/Cell.h"
#include "redblood/Interpolation.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace redblood
  {
    namespace
    {
      template<class T>
      int figureNearness(DivideConquer<T> &dnc, LatticeVector const &key,
                         LatticePosition const &vertex, PhysicalDistance const &haloLength)
      {
        if (haloLength + haloLength > dnc.GetBoxSize())
        {
          return 0;
        }

        int result = 0;

        for (size_t d(1); d < (1 << 6); d <<= 1)
        {
          LatticePosition const translated(CellReference::directions(d) * haloLength);

          if (not (key == dnc.DowngradeKey(vertex + translated)))
          {
            result |= d;
          }
        }

        return result;
      }

      template<class T>
      CellReference initCellRef(DivideConquer<T> &dnc, CellContainer::const_iterator cellid,
                                site_t nodeid, LatticeVector const &key,
                                LatticePosition const &vertex, PhysicalDistance const &haloLength)
      {
        int const isNearBorder = figureNearness(dnc, key, vertex, haloLength);
        CellReference result = { cellid, nodeid, isNearBorder };
        return result;
      }

      void initializeCells(DivideConquer<CellReference> &dnc, MeshData::Vertices const &vertices,
                           CellContainer::const_iterator cellid, PhysicalDistance haloLength)
      {
        typedef DivideConquer<CellReference> DnC;
        typedef DnC::key_type key_type;
        typedef MeshData::Vertices::const_iterator vertex_iterator;
        vertex_iterator i_first = vertices.begin();
        vertex_iterator const i_end = vertices.end();

        for (site_t i(0); i_first != i_end; ++i_first, ++i)
        {
          key_type const key = dnc.DowngradeKey(*i_first);
          dnc.insert(key, initCellRef(dnc, cellid, i, key, *i_first, haloLength));
        }
      }

      // avoids a warning
#     ifndef HEMELB_DOING_UNITTESTS
      void initializeCells(DivideConquer<CellReference> &dnc, CellContainer const &cells,
                           PhysicalDistance haloLength)
      {
        CellContainer::const_iterator i_first = cells.begin();
        CellContainer::const_iterator const i_end = cells.end();

        for (site_t i(0); i_first != i_end; ++i_first, ++i)
        {
          initializeCells(dnc, (*i_first)->GetVertices(), i_first, haloLength);
        }
      }
#     endif

      // Compare distance between vertices
      template<class T_FUNCTION>
      bool nextDistance(T_FUNCTION const &strictlyLarger, DivideConquerCells::const_iterator &first,
                        DivideConquerCells::const_iterator const &end,
                        DivideConquerCells::const_iterator const &main, PhysicalDistance dist)
      {
        auto const mainCell = main.GetCell();
        typedef DivideConquerCells::const_iterator cit;
        auto goodCellPair = [&mainCell, &strictlyLarger](cit const &i)
        {
          return strictlyLarger(i.GetCell(), mainCell);
        };
        auto goodDistance = [&main, &dist](cit const &i)
        {
          return (*main - *i).GetMagnitude() < dist;
        };
        for (; first != end; ++first)
        {
          if (goodCellPair(first) and goodDistance(first))
          {
            return true;
          }
        }
        return false;
      }

      // avoids a warning
#     ifndef HEMELB_DOING_UNITTESTS
      void spreadForce(LatticePosition const &node, geometry::LatticeData &latticeData,
                       stencil::types stencil, LatticeForceVector const &force)
      {
        proc_t procid;
        site_t siteid;
        InterpolationIterator spreader = interpolationIterator(node, stencil);

        for (; spreader; ++spreader)
        {
          if (latticeData.GetContiguousSiteId(*spreader, procid, siteid))
          {
            latticeData.GetSite(siteid).AddToForce(force * spreader.weight());
          }
        }
      }
#     endif
    } // anonymous namespace
#ifndef HEMELB_DOING_UNITTESTS
    //! Constructor
    DivideConquerCells::DivideConquerCells(CellContainer const &cells, PhysicalDistance boxsize,
                                           PhysicalDistance halosize) :
        DivideConquer<CellReference>(boxsize), haloLength(halosize), cells(cells)
    {
      initializeCells(*static_cast<base_type *>(this), GetCells(), haloLength);
    }

    void DivideConquerCells::update()
    {
      iterator i_first = begin();
      iterator const i_end = end();

      for (; i_first != i_end;)
      {
        key_type const key = base_type::DowngradeKey(*i_first);
        i_first.GetCellReference().isNearBorder = figureNearness(*this, key, *i_first, haloLength);

        if (not (key == i_first.GetKey()))
        {
          base_type::insert(key, i_first.GetCellReference());
          base_type::iterator const next_iter = base_type::erase((base_type::iterator) i_first);
          i_first = iterator(*this, next_iter);
        }
        else
        {
          ++i_first;
        }
      }
    }

    DivideConquerCells::const_range DivideConquerCells::operator()(LatticeVector const &pos) const
    {
      base_type::const_range const boxrange = base_type::equal_range(pos);
      return const_range(const_iterator(*this, boxrange.first),
                         const_iterator(*this, boxrange.second));
    }

    bool DivideConquerCells::pair_range::nextDist()
    {
      typedef decltype(owner.cells) Cells;
      typedef Cells::const_reference Input;
      auto strictly_less = Cells::key_compare();
      auto strictlyLarger = [&strictly_less](Input _a, Input _b)
      {
        return _a != _b and not strictly_less(_a, _b);
      };
      return nextDistance(strictlyLarger, currents.second, ends.second, currents.first, maxdist);
    }

    bool DivideConquerCells::pair_range::doBox()
    {
      LatticeVector const key(box == CellReference::NONE ?
        currents.first.GetKey() :
        currents.first.GetKey() + CellReference::idirections(box));
      DivideConquerCells::const_range const boxits = owner(key);

      if (box == CellReference::NONE)
      {
        currents.second = currents.first;
        ++currents.second;
      }
      else
      {
        currents.second = boxits.first;
      }

      ends.second = boxits.second;
      return nextDist();
    }

    bool DivideConquerCells::pair_range::operator++()
    {
      if (not is_valid())
      {
        return false;
      }

      // First try and finds next pair in current range
      if (currents.second != ends.second)
      {
        ++currents.second;

        if (nextDist())
        {
          return true;
        }
      }

      // If reaches here, then should check which box we are currently doing
      if (currents.first.GetNearBorder())
      {
        if (box)
        {
          box = CellReference::Borders(int(box) << 1);
        }
        else
        {
          box = CellReference::Borders(1);
        }

        while (box < CellReference::LAST)
        {
          if (doBox())
          {
            return true;
          }

          box = CellReference::Borders(int(box) << 1);
        }
      }

      // If reaches here, then should increment main iterator and start with same
      // box
      if (++currents.first == ends.first)
      {
        return false;
      }

      box = CellReference::NONE;
      return doBox() ?
        true :
        operator++();
    }

    DivideConquerCells::pair_range::pair_range(DivideConquerCells const &owner,
                                               iterator const &begin, iterator const &end,
                                               PhysicalDistance maxdist) :
        maxdist(maxdist), box(CellReference::NONE), currents(begin, end), ends(end, end),
            owner(owner)
    {
      // No throw garantee. Makes iterator invalid instead.
      try
      {
        // Could be invalid from start
        if (not is_valid())
        {
          return;
        }

        // Iterates to first valid item, if any
        if (not doBox())
        {
          operator++();
        }
      }
      catch (std::exception const &e)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("*** Encountered error while initializing pair iterator: %s\n",
                                                      e.what());
        currents.first = ends.first;
      }
      catch (...)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("*** Encountered error while initializing pair iterator.");
        currents.first = ends.first;
      }
    }

    DivideConquerCells::pair_range DivideConquerCells::pair_begin(PhysicalDistance maxdist) const
    {
      return pair_range(*this, begin(), end(), maxdist);
    }

    //! Computes cell <-> cell interactions and spread to grid
    void addCell2CellInteractions(DivideConquerCells const &dnc, Node2NodeForce const &functional,
                                  stencil::types stencil, geometry::LatticeData &latticeData)
    {
      DivideConquerCells::pair_range range(dnc.pair_begin(functional.cutoff));

      for (; range.is_valid(); ++range)
      {
        LatticeForceVector const force(functional(*range->first, *range->second));
        // spread to the grid from from one node and from the other
        spreadForce(*range->first, latticeData, stencil, force);
        spreadForce(*range->second, latticeData, stencil, -force);
      }
    }

    bool DivideConquerCells::insert(CellContainer::value_type cell)
    {
      auto inserted = cells.insert(cell);
      if (inserted.second)
        initializeCells(*this, (*inserted.first)->GetVertices(), inserted.first, haloLength);
      return true;
    }

    bool DivideConquerCells::remove(CellContainer::value_type cell)
    {
      auto const cellIterator = cells.find(cell);
      if (cellIterator == cells.end())
        return false;
      auto i_first = base_type::cbegin();
      auto const i_end = base_type::cend();
      while (i_first != i_end)
      {
        auto const i_current = i_first;
        // Move to next element *before* possibly erasing current element
        // Not sure whether i_first can be incremented after an erase.
        ++i_first;
        if (*i_current->second.cellIterator == cell)
        {
          base_type::erase(i_current);
        }
      }
      cells.erase(cellIterator);
      return true;
    }
#endif
  }
} // namespace hemelb::redblood
