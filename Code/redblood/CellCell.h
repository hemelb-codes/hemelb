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

#include "units.h"
#include "Exception.h"
#include "geometry/LatticeData.h"
#include "redblood/DivideConquer.h"
#include "redblood/Cell.h"
#include "redblood/Node2Node.h"
#include "redblood/stencil.h"

namespace hemelb
{
  namespace redblood
  {
    struct CellReference
    {
      //! Index of cell in input container
      site_t cellIndex;
      //! Index of node in mesh
      site_t nodeIndex;
      //! Whether the node is near the border of the cube;
      int isNearBorder;

      //! Names for each border
      enum Borders
      {
        NONE = 0,
        TOP = 1,
        BOTTOM = 2,
        NORTH = 4,
        SOUTH = 8,
        WEST = 16,
        EAST = 32,
        LAST = 64
      };

      static LatticeVector idirections(Borders border)
      {
        switch (border)
        {
          case TOP:
            return LatticeVector(1, 0, 0);

          case BOTTOM:
            return LatticeVector(-1, 0, 0);

          case NORTH:
            return LatticeVector(0, 1, 0);

          case SOUTH:
            return LatticeVector(0, -1, 0);

          case WEST:
            return LatticeVector(0, 0, -1);

          case EAST:
            return LatticeVector(0, 0, 1);

          default:
            return LatticeVector(0, 0, 0);
        };
      }
      static LatticePosition directions(Borders border)
      {
        return LatticePosition(idirections(border));
      }

      static LatticeVector idirections(size_t border)
      {
        assert(border == TOP or border == BOTTOM or border == NORTH or border == SOUTH
               or border == EAST or border == WEST);
        return idirections(Borders(border));
      }
      static LatticePosition directions(size_t border)
      {
        return LatticePosition(idirections(border));
      }
    };

    //! Organizes nodes in cells in boxes
    //! The object is to easily check nodes that are within interaction distance
    class DivideConquerCells : protected DivideConquer<CellReference>
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
      DivideConquerCells(CellContainer const &cells, PhysicalDistance boxsize,
                         PhysicalDistance halosize);

      //! Gets all nodes in a box
      const_range operator()(LatticeVector const &pos) const;
      //! Gets all nodes in a box
      const_range operator()(LatticePosition const &pos) const
      {
        return this->operator()(base_type::DowngradeKey(pos));
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
      size_t size() const
      {
        return base_type::size();
      }

      //! Loops over pair of vertices closer than input distance
      pair_range pair_begin(PhysicalDistance maxdist) const;

      //! Distance from border below which an object is in the halo
      PhysicalDistance GetHaloLength() const
      {
        return haloLength;
      }
      //! Size of each box
      PhysicalDistance GetBoxSize() const
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

      protected:
      //! Distance from border below which an object is in the halo
      PhysicalDistance const haloLength;
      //! Container of cells
      CellContainer cells;
    };

    class DivideConquerCells::pair_range
    {
      //! Parent iterator
      typedef DivideConquerCells::const_iterator iterator;

      public:
      typedef std::pair<iterator, iterator> value_type;

      //! Constructs a pair range iterator
      pair_range(DivideConquerCells const &owner, iterator const &begin, iterator const &end,
                 PhysicalDistance maxdist);

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
        return &(this->operator*());
      }

      //! Goes to next value and returns true if valid
      bool operator++();

      protected:
      //! Maximum distance for which to report pair
      PhysicalDistance maxdist;
      //! Current box we are working on
      CellReference::Borders box;
      //! Iterator for main item
      value_type currents;
      //! range for iteration over second item
      value_type ends;
      //! Container over which this one iterates
      DivideConquerCells const &owner;

      bool doBox();
      bool nextDist();
    };

    //! Computes cell <-> cell interactions and spread to grid
    //! Given a partition of the cells' nodes and node <-> node interaction
    //! functional, computes the short-range that can occur between cells that are
    //! too close to one another. The interaction forces are computed and spread to
    //! the lattice.
    void addCell2CellInteractions(DivideConquerCells const &dnc, Node2NodeForce const &functional,
                                  stencil::types stencil, geometry::LatticeData &latticeData);
  }
}  // hemelb::redblood

#endif
