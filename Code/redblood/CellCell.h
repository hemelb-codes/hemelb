//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_CELL_CELL_INTERACTION_H
#define HEMELB_REDBLOOD_CELL_CELL_INTERACTION_H

#include <vector>
#include <assert.h>
#include <memory>
#include <initializer_list>
#include <type_traits>

#include "units.h"
#include "Exception.h"
#include "geometry/LatticeData.h"
#include "redblood/DivideConquer.h"
#include "redblood/Cell.h"
#include "redblood/Node2Node.h"
#include "redblood/stencil.h"
#include "redblood/Interpolation.h"
#include "redblood/Borders.h"
#include "redblood/parallel/CellParallelization.h"

namespace hemelb
{
  namespace redblood
  {
    // References a node of a mesh in the divide-and-conquer box
    class CellReference
    {
      public:
        //! Index of cell in input container
        CellContainer::const_iterator cellIterator;
        //! Index of node in mesh
        site_t nodeIndex;
        //! Id of the nearest borders
        size_t nearBorder;
    };

    //! Organizes nodes in cells in boxes
    //! The object is to easily check nodes that are within interaction distance
    class DivideConquerCells : public DivideConquer<CellReference>
    {
        //! Type of the base class
        typedef DivideConquer<CellReference> base_type;

      public:
        //! Iterates over vertices
        //! Wraps a multimap iterator. As such it provides the exact same
        //! garantees.
        class iterator;
        //! Iterates over vertices
        //! Wraps a multimap iterator. As such it provides the exact same
        //! garantees.
        class const_iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        //! Iterator pair over single box
        typedef std::pair<const_iterator, const_iterator> const_range;
        //! Iterates over pairs of nodes that are close to one another
        //! Each pair is visited only once.
        class pair_range;

        //! Constructor
        DivideConquerCells(CellContainer const &cells, LatticeDistance boxsize,
                           LatticeDistance halosize);

        //! Gets all nodes in a box
        const_range operator()(LatticeVector const &pos) const;
        //! Gets all nodes in a box
        const_range operator()(LatticePosition const &pos) const
        {
          return this->operator()(base_type::DowngradeKey(pos));
        }
        template<class T>
        typename std::enable_if<std::is_integral<T>::value, const_range>::type operator()(
            std::initializer_list<T> pos) const
        {
          assert(pos.size() == 3);
          return operator()(LatticeVector(*pos.begin(),
                                          *std::next(pos.begin()),
                                          *std::next(std::next(pos.begin()))));
        }
        template<class T>
        typename std::enable_if<std::is_floating_point<T>::value, const_range>::type operator()(
            std::initializer_list<T> pos) const
        {
          assert(pos.size() == 3);
          return operator()(LatticePosition(*pos.begin(),
                                            *std::next(pos.begin()),
                                            *std::next(std::next(pos.begin()))));
        }

// Implementation of DivideConquerCell::iterator
#define HEMELB_DOING_NONCONST
#include "CellCellIterator.impl.h"
#undef HEMELB_DOING_NONCONST
// Implementation of DivideConquerCell::const_iterator
#include "CellCellIterator.impl.h"

        iterator begin()
        {
          return iterator(*this, base_type::begin());
        }
        iterator end()
        {
          return iterator(*this, base_type::end());
        }
        const_iterator begin() const
        {
          return const_iterator(*this, base_type::begin());
        }
        const_iterator end() const
        {
          return const_iterator(*this, base_type::end());
        }
        reverse_iterator rbegin()
        {
          return reverse_iterator(end());
        }
        reverse_iterator rend()
        {
          return reverse_iterator(begin());
        }
        const_reverse_iterator rbegin() const
        {
          return const_reverse_iterator(end());
        }
        const_reverse_iterator rend() const
        {
          return const_reverse_iterator(begin());
        }
        base_type::size_type size() const
        {
          return base_type::size();
        }

        //! Loops over pair of vertices closer than input distance
        pair_range pair_begin(LatticeDistance maxdist) const;

        //! Distance from border below which an object is in the halo
        LatticeDistance GetHaloLength() const
        {
          return haloLength;
        }
        //! Size of each box
        LatticeDistance GetBoxSize() const
        {
          return base_type::GetBoxSize();
        }
        //! Cells that are present in this object
        CellContainer const &GetCells() const
        {
          return cells;
        }

        //! After vertices have moved, update mapping and whether it is near
        //! boundary
        void update();
        //! recomputes using current cells
        void SetBoxSizeAndHalo(LatticeDistance boxSize, LatticeDistance halo);

        //! In a parallel simulation, some cells will leave/enter the domain.
        //! \param[in] 3-tuple with the newly owned cells, the disowned cells, and the lent cells
        void update(parallel::ExchangeCells::ChangedCells const& changedCells);

        //! Insert a new cell
        //! Returns true if the cell was inserted, false if it already existed.
        bool insert(CellContainer::value_type cell);
        //! Removes a cell
        //! Returns true if the cell did exist.
        bool remove(CellContainer::value_type cell);

      protected:
        //! Distance from border below which an object is in the halo
        LatticeDistance haloLength;
        //! Container of cells
        CellContainer cells;
        //! Keeps track of which cells in the box are lent (currentlyLentCells must always be a subset of DivideConquerCells::cells)
        CellContainer currentlyLentCells;
        //! Helper method that combines all the CellContainers in a LentCells object into a single CellContainer
        CellContainer LentCellsToSingleContainer(parallel::CellParallelization::LentCells const &lentCells) const;
    };

    class DivideConquerCells::pair_range
    {
        //! Parent iterator
        typedef DivideConquerCells::const_iterator iterator;

      public:
        typedef std::pair<iterator, iterator> value_type;

        //! Constructs a pair range iterator
        pair_range(DivideConquerCells const &owner, iterator const &begin, iterator const &end,
                   LatticeDistance maxdist);

        //! Whether current iteration is valid
        bool is_valid() const
        {
          return currents.first != ends.first;
        }
        //! Whether current iteration is valid
        operator bool() const
        {
          return is_valid();
        }

        value_type const &operator*() const
        {
#ifndef NDEBUG

          if (not is_valid())
          {
            throw Exception() << "Iterator is invalid\n";
          }

#endif
          return currents;
        }
        value_type const *operator->() const
        {
          return & (this->operator*());
        }

        //! Goes to next value and returns true if valid
        bool operator++();

      protected:
        //! Maximum distance for which to report pair
        LatticeDistance maxdist;
        //! Iterates over boxes we want to see
        BorderBoxIterator box_iterator;
        //! Iterator for main item
        value_type currents;
        //! range for iteration over second item
        value_type ends;
        //! Container over which this one iterates
        DivideConquerCells const &owner;

        bool doBox();
        bool nextDist();
    };

    namespace
    {
      //! Spread force from given vertices to lattice-sites
      template<class STENCIL>
      void spreadForce(LatticePosition const &vertex, geometry::LatticeData &latticeData,
                       LatticeForceVector const &force)
      {
        proc_t procid;
        site_t siteid;
        InterpolationIterator<STENCIL> spreader = interpolationIterator<STENCIL>(vertex);

        for (; spreader; ++spreader)
        {
          if (latticeData.GetContiguousSiteId(*spreader, procid, siteid))
          {
            latticeData.GetSite(siteid).AddToForce(force * spreader.weight());
          }
        }
      }
    }
    //! \brief Computes cell <-> cell interactions and spread to grid
    //! \details Given a partition of the cells' nodes and node <-> node interaction
    //! functional, computes the short-range that can occur between cells that are
    //! too close to one another. The interaction forces are computed and spread to
    //! the lattice.
    template<class STENCIL>
    void addCell2CellInteractions(DivideConquerCells const &dnc, Node2NodeForce const &functional,
                                  geometry::LatticeData &latticeData)
    {
      auto range = dnc.pair_begin(functional.cutoff);

      for (; range.is_valid(); ++range)
      {
        LatticeForceVector const force(functional(*range->first, *range->second));
        // spread to the grid from from one node and from the other
        spreadForce<STENCIL>(*range->first, latticeData, force);
        spreadForce<STENCIL>(*range->second, latticeData, -force);
      }
    }

  }
} // hemelb::redblood

#endif
