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

namespace hemelb { namespace redblood {

namespace {
  template<class T> int figure_nearness(
      DivideConquer<T> &_dnc,
      LatticeVector const &_key,
      LatticePosition const &_vertex,
      PhysicalDistance const &_haloLength
  ) {
    if(_haloLength + _haloLength > _dnc.GetBoxSize()) return 0;
    int result = 0;
    for(size_t d(1); d < (1 << 6); d <<= 1) {
      LatticePosition const translated(
         CellReference::directions(d) * _haloLength
      );
      if(not (_key == _dnc.DowngradeKey(_vertex + translated)))
        result |= d;
    }
    return result;
  }

  template<class T> CellReference init_cell_ref(
      DivideConquer<T> &_dnc,
      site_t _cellid, site_t _nodeid,
      LatticeVector const &_key,
      LatticePosition const &_vertex,
      PhysicalDistance const &_haloLength
  ) {
    int const isNearBorder = figure_nearness(_dnc, _key, _vertex, _haloLength);
    CellReference result = {_cellid, _nodeid, isNearBorder};
    return result;
  }

  void initialize_cells(
      DivideConquer<CellReference> &_dnc,
      MeshData::Vertices const &_vertices,
      site_t _cellid,
      PhysicalDistance _haloLength ) {
    typedef DivideConquer<CellReference> DnC;
    typedef DnC::key_type key_type;
    typedef MeshData::Vertices::const_iterator vertex_iterator;
    vertex_iterator i_first = _vertices.begin();
    vertex_iterator const i_end = _vertices.end();
    for(site_t i(0); i_first != i_end; ++i_first, ++i) {
      key_type const key = _dnc.DowngradeKey(*i_first);
      _dnc.insert(
          key,
          init_cell_ref(_dnc, _cellid, i, key, *i_first, _haloLength)
      );
    }
  }

  void initialize_cells(
      DivideConquer<CellReference> &_dnc,
      CellContainer const& _cells,
      PhysicalDistance _haloLength ) {
    CellContainer::const_iterator i_first = _cells.begin();
    CellContainer::const_iterator const i_end = _cells.end();
    for(site_t i(0); i_first != i_end; ++i_first, ++i)
      initialize_cells(_dnc, (*i_first)->GetVertices(), i, _haloLength);
  }

  // Compare distance between vertices
  bool next_dist(
      DivideConquerCells::const_iterator &_first,
      DivideConquerCells::const_iterator const & _end,
      DivideConquerCells::const_iterator const & _main,
      PhysicalDistance _dist
  ) {
    for(; _first != _end; ++_first)
      if(_first.GetCellIndex() > _main.GetCellIndex()
          and (*_main - *_first).GetMagnitude() < _dist)
        return true;
    return false;
  }

  void spreadForce(
      LatticePosition const &_node,
      geometry::LatticeData &_latticeData,
      stencil::types _stencil,
      LatticeForceVector const &_force
  ) {
    proc_t procid;
    site_t siteid;
    InterpolationIterator spreader = interpolationIterator(_node, _stencil);
    for(; spreader; ++spreader)
      if(_latticeData.GetContiguousSiteId(*spreader, procid, siteid))
        _latticeData.GetSite(siteid).AddToForce(_force * spreader.weight());
  }
}

#ifndef HEMELB_DOING_UNITTESTS
//! Constructor
DivideConquerCells :: DivideConquerCells(
    CellContainer const &_cells,
    PhysicalDistance _boxsize, PhysicalDistance _halosize
) : DivideConquer<CellReference>(_boxsize),
    haloLength_(_halosize), cells_(_cells) {
  try {
    initialize_cells(*static_cast<base_type*>(this), _cells, haloLength_);
  } catch(...) {}
}

void DivideConquerCells :: update() {
   iterator i_first = begin();
   iterator const i_end = end();
   for(; i_first != i_end; ++i_first) {
      key_type const key = base_type::DowngradeKey(*i_first);
      i_first.GetCellReference().isNearBorder = figure_nearness(
          *this, key, *i_first, haloLength_
      );
      if(not (key == i_first.GetKey())) {
        base_type::insert(key, i_first.GetCellReference());
        base_type::erase((base_type::iterator) i_first);
      }
   }
}

DivideConquerCells::const_range DivideConquerCells::operator()(
    LatticeVector const &_pos) const {
  base_type::const_range const boxrange = base_type::equal_range(_pos);
  return const_range(
      const_iterator(*this, boxrange.first),
      const_iterator(*this, boxrange.second)
  );
}

bool DivideConquerCells::pair_range::next_dist_() {
  return next_dist(currents_.second, ends_.second, currents_.first, maxdist_);
}

bool DivideConquerCells::pair_range::do_box_() {
  LatticeVector const key(
      box_ == CellReference::NONE ?
        currents_.first.GetKey():
        currents_.first.GetKey() + CellReference::idirections(box_)
  );
  DivideConquerCells::const_range const boxits = owner_(key);
  if(box_ == CellReference::NONE) {
    currents_.second = currents_.first;
    ++currents_.second;
  } else
    currents_.second = boxits.first;
  ends_.second = boxits.second;
  return next_dist_();
}

bool DivideConquerCells::pair_range::operator++() {

  if(not is_valid()) return false;

  // First try and finds next pair in current range
  if(currents_.second != ends_.second) {
    ++currents_.second;
    if(next_dist_()) return true;
  }

  // If reaches here, then should check which box we are currently doing
  if(currents_.first.GetNearBorder()) {
    if(box_) box_ = CellReference::Borders(int(box_) << 1);
    else box_ = CellReference::Borders(1);
    while(box_ < CellReference::LAST) {
      if(do_box_()) return true;
      box_ = CellReference::Borders(int(box_) << 1);
    }
  }

  // If reaches here, then should increment main iterator and start with same
  // box
  if(++currents_.first == ends_.first) return false;
  box_ = CellReference::NONE;
  return do_box_() ? true: operator++();
}

DivideConquerCells::pair_range::pair_range(
    DivideConquerCells const &_owner,
    iterator const &_begin,
    iterator const &_end,
    PhysicalDistance _maxdist
) : maxdist_(_maxdist), box_(CellReference::NONE),
    currents_(_begin, _end), ends_(_end, _end), owner_(_owner) {
  // No throw garantee. Makes iterator invalid instead.
  try {
    // Could be invalid from start
    if(not is_valid()) return;
    // Iterates to first valid item, if any
    if(not do_box_()) operator++();
  } catch(std::exception const &_e) {
    log::Logger::Log<log::Debug, log::OnePerCore>(
      "*** Encountered error while initializing pair iterator: %s\n",
      _e.what()
    );
    currents_.first = ends_.first;
  } catch(...) {
    log::Logger::Log<log::Debug, log::OnePerCore>(
      "*** Encountered error while initializing pair iterator.");
    currents_.first = ends_.first;
  }
}

DivideConquerCells::pair_range DivideConquerCells::pair_begin(
    PhysicalDistance _maxdist) const {
  return pair_range(*this, begin(), end(), _maxdist);
}

//! Computes cell <-> cell interactions and spread to grid
void addCell2CellInteractions(
    DivideConquerCells const &_dnc,
    Node2NodeForce const &_functional,
    stencil::types _stencil,
    geometry::LatticeData &_latticeData
) {
  DivideConquerCells :: pair_range range(_dnc.pair_begin(_functional.cutoff));
  for(; range.is_valid(); ++range) {
    LatticeForceVector const force(_functional(*range->first, *range->second));
    // spread to the grid from from one node and from the other
    spreadForce(*range->first, _latticeData, _stencil, force);
    spreadForce(*range->second, _latticeData, _stencil, -force);
  }
}
#endif

}} // namespace hemelb::redblood
