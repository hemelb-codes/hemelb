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
      DivideConquer<T> &dnc,
      LatticeVector const &key,
      LatticePosition const &vertex,
      PhysicalDistance const &haloLength
  ) {
    if(haloLength + haloLength > dnc.GetBoxSize()) return 0;
    int result = 0;
    for(size_t d(1); d < (1 << 6); d <<= 1) {
      LatticePosition const translated(
         CellReference::directions(d) * haloLength
      );
      if(not (key == dnc.DowngradeKey(vertex + translated)))
        result |= d;
    }
    return result;
  }

  template<class T> CellReference init_cell_ref(
      DivideConquer<T> &dnc,
      site_t cellid, site_t nodeid,
      LatticeVector const &key,
      LatticePosition const &vertex,
      PhysicalDistance const &haloLength
  ) {
    int const isNearBorder = figure_nearness(dnc, key, vertex, haloLength);
    CellReference result = {cellid, nodeid, isNearBorder};
    return result;
  }

  void initialize_cells(
      DivideConquer<CellReference> &dnc,
      MeshData::Vertices const &vertices,
      site_t cellid,
      PhysicalDistance haloLength ) {
    typedef DivideConquer<CellReference> DnC;
    typedef DnC::key_type key_type;
    typedef MeshData::Vertices::const_iterator vertex_iterator;
    vertex_iterator i_first = vertices.begin();
    vertex_iterator const i_end = vertices.end();
    for(site_t i(0); i_first != i_end; ++i_first, ++i) {
      key_type const key = dnc.DowngradeKey(*i_first);
      dnc.insert(
          key,
          init_cell_ref(dnc, cellid, i, key, *i_first, haloLength)
      );
    }
  }

  void initialize_cells(
      DivideConquer<CellReference> &dnc,
      CellContainer const& cells,
      PhysicalDistance haloLength ) {
    CellContainer::const_iterator i_first = cells.begin();
    CellContainer::const_iterator const i_end = cells.end();
    for(site_t i(0); i_first != i_end; ++i_first, ++i)
      initialize_cells(dnc, (*i_first)->GetVertices(), i, haloLength);
  }

  // Compare distance between vertices
  bool next_dist(
      DivideConquerCells::const_iterator &first,
      DivideConquerCells::const_iterator const & end,
      DivideConquerCells::const_iterator const & main,
      PhysicalDistance dist
  ) {
    for(; first != end; ++first)
      if(first.GetCellIndex() > main.GetCellIndex()
          and (*main - *first).GetMagnitude() < dist)
        return true;
    return false;
  }

  void spreadForce(
      LatticePosition const &node,
      geometry::LatticeData &latticeData,
      stencil::types stencil,
      LatticeForceVector const &force
  ) {
    proc_t procid;
    site_t siteid;
    InterpolationIterator spreader = interpolationIterator(node, stencil);
    for(; spreader; ++spreader)
      if(latticeData.GetContiguousSiteId(*spreader, procid, siteid))
        latticeData.GetSite(siteid).AddToForce(force * spreader.weight());
  }
}

#ifndef HEMELB_DOING_UNITTESTS
//! Constructor
DivideConquerCells :: DivideConquerCells(
    CellContainer const &cells,
    PhysicalDistance boxsize, PhysicalDistance halosize
) : DivideConquer<CellReference>(boxsize),
    haloLength_(halosize), cells_(cells) {
  try {
    initialize_cells(*static_cast<base_type*>(this), cells, haloLength_);
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
    LatticeVector const &pos) const {
  base_type::const_range const boxrange = base_type::equal_range(pos);
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
    DivideConquerCells const &owner,
    iterator const &begin,
    iterator const &end,
    PhysicalDistance maxdist
) : maxdist_(maxdist), box_(CellReference::NONE),
    currents_(begin, end), ends_(end, end), owner_(owner) {
  // No throw garantee. Makes iterator invalid instead.
  try {
    // Could be invalid from start
    if(not is_valid()) return;
    // Iterates to first valid item, if any
    if(not do_box_()) operator++();
  } catch(std::exception const &e) {
    log::Logger::Log<log::Debug, log::OnePerCore>(
      "*** Encountered error while initializing pair iterator: %s\n",
      e.what()
    );
    currents_.first = ends_.first;
  } catch(...) {
    log::Logger::Log<log::Debug, log::OnePerCore>(
      "*** Encountered error while initializing pair iterator.");
    currents_.first = ends_.first;
  }
}

DivideConquerCells::pair_range DivideConquerCells::pair_begin(
    PhysicalDistance maxdist) const {
  return pair_range(*this, begin(), end(), maxdist);
}

//! Computes cell <-> cell interactions and spread to grid
void addCell2CellInteractions(
    DivideConquerCells const &dnc,
    Node2NodeForce const &functional,
    stencil::types stencil,
    geometry::LatticeData &latticeData
) {
  DivideConquerCells :: pair_range range(dnc.pair_begin(functional.cutoff));
  for(; range.is_valid(); ++range) {
    LatticeForceVector const force(functional(*range->first, *range->second));
    // spread to the grid from from one node and from the other
    spreadForce(*range->first, latticeData, stencil, force);
    spreadForce(*range->second, latticeData, stencil, -force);
  }
}
#endif

}} // namespace hemelb::redblood
