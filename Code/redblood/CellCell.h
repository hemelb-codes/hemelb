//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_CELL_CELL_INTERACTION_H_
#define HEMELB_REDBLOOD_CELL_CELL_INTERACTION_H_

#include <vector>
#include "units.h"
#include "redblood/Cell.h"
#include "redblood/DivideConquer.h"


namespace hemelb { namespace redblood {

struct CellReference {
  //! Index of cell in input container
  site_t cellIndex;
  //! Index of node in mesh
  site_t nodeIndex;
  //! Whether the node is near the border of the cube;
  int isNearBorder;

  //! Names for each border
  enum Borders {
    TOP = 1,
    BOTTOM = 2,
    NORTH = 4,
    SOUTH = 8,
    WEST = 16,
    EAST = 32
  };

  static LatticePosition directions(Borders _border) {
    switch(_border) {
      case TOP:    return LatticePosition(1, 0, 0);
      case BOTTOM: return LatticePosition(-1, 0, 0);
      case NORTH:   return LatticePosition(0, 1, 0);
      case SOUTH:  return LatticePosition(0, -1, 0);
      case WEST:  return LatticePosition(0, 0, -1);
      case EAST:   return LatticePosition(0, 0, 1);
    };
  }
  static LatticePosition directions(size_t _border) {
    assert(
        _border == TOP or _border == BOTTOM or _border == NORTH
        or _border == SOUTH or _border == EAST or  _border == WEST
    );
    return directions(Borders(_border));
  }
};

//! Contains all cells...
typedef std::vector<Cell> CellContainer;

//! Organizes nodes in cells in boxes
//! The object is to easily check nodes that are within interaction distance
class DivideConquerCells : protected DivideConquer<CellReference> {
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

    //! Constructor
    DivideConquerCells(
        boost::shared_ptr<CellContainer> const &_cells,
        PhysicalDistance _boxsize, PhysicalDistance _halosize
    );

    //! Gets all nodes in a box
    const_range operator()(LatticeVector const &_pos) const;
    //! Gets all nodes in a box
    const_range operator()(LatticePosition const &_pos) const {
      return this->operator()(base_type::DowngradeKey(_pos));
    }

    // Implementation of DivideConquerCell::iterator
#   define HEMELB_DOING_NONCONST
#     include "CellCellIterator.impl.h"
#   undef HEMELB_DOING_NONCONST
    // Implementation of DivideConquerCell::const_iterator
#   include "CellCellIterator.impl.h"

    iterator begin() { return iterator(*this, base_type::begin()); }
    iterator end() { return iterator(*this, base_type::end()); }
    const_iterator begin() const {
      return const_iterator(*this, base_type::begin());
    }
    const_iterator end() const {
      return const_iterator(*this, base_type::end());
    }
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rbegin() const {
      return const_reverse_iterator(end());
    }
    const_reverse_iterator rend() const {
      return const_reverse_iterator(begin());
    }
    size_t size() const { return base_type::size(); }

    //! Distance from border below which an object is in the halo
    PhysicalDistance GetHaloLength() const { return haloLength_; }
    //! Size of each box
    PhysicalDistance GetBoxSize() const { return base_type::GetBoxSize(); }
    //! Cells that are present in this object
    boost::shared_ptr<CellContainer const> GetCells() const { return cells_; }

    //! After vertices have moved, update mapping and whether it is near
    //! boundary
    void update();
  protected:
    //! Distance from border below which an object is in the halo
    PhysicalDistance const haloLength_;
    //! Container of cells
    boost::shared_ptr<CellContainer> cells_;
};


}} // hemelb::redblood

#endif
