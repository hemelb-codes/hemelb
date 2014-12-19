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
      MeshData::t_Vertices const &_vertices,
      site_t _cellid,
      PhysicalDistance _haloLength ) {
    typedef DivideConquer<CellReference> DnC;
    typedef DnC::key_type key_type;
    typedef MeshData::t_Vertices::const_iterator vertex_iterator;
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
      boost::shared_ptr<CellContainer const> _cells,
      PhysicalDistance _haloLength ) {
    CellContainer::const_iterator i_first = _cells->begin();
    CellContainer::const_iterator const i_end = _cells->end();
    for(site_t i(0); i_first != i_end; ++i_first, ++i)
      initialize_cells(_dnc, i_first->GetVertices(), i, _haloLength);
  }
}

#ifndef HEMELB_DOING_UNITTESTS
//! Constructor
DivideConquerCells :: DivideConquerCells(
    boost::shared_ptr<CellContainer> const &_cells,
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
#endif

}} // namespace hemelb::redblood
