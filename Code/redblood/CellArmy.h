//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELL_ARMY_H
#define HEMELB_UNITTESTS_REDBLOOD_CELL_ARMY_H

#include <algorithm>
#include <vector>
#include <memory>

#include "redblood/Cell.h"
#include "redblood/CellCell.h"
#include "geometry/LatticeData.h"

namespace hemelb { namespace redblood {

//! \brief Federates the cells together so we can apply ops simultaneously
template<class KERNEL> class CellArmy {
  public:
    //! Interaction terms between cells
    Node2NodeForce cell2Cell;
    //! Stencil
    stencil::types stencil = stencil::FOUR_POINT;

    CellArmy(geometry::LatticeData &_latDat,
        CellContainer const &_cells,
        PhysicalDistance _boxsize=10.0, PhysicalDistance _halo=2.0)
      : latticeData_(_latDat), cells_(_cells),
      dnc_(_cells, _boxsize, _halo) {}

    //! Performs fluid to lattice interactions
    void fluid2CellInteractions();

    //! Performs lattice to fluid interactions
    void cell2FluidInteractions();

#   ifdef HEMELB_DOING_UNITTESTS
      //! Updates divide and conquer
      void updateDNC() { dnc_.update(); }
#   endif

  protected:
    //! All lattice information and then some
    geometry::LatticeData &latticeData_;
    //! Contains all cells
    CellContainer cells_;
    //! Divide and conquer object
    DivideConquerCells dnc_;
    //! A work array with forces/positions
    std::vector<LatticePosition> work_;
};

template<class KERNEL>
  void CellArmy<KERNEL> :: fluid2CellInteractions() {
    std::vector<LatticePosition> & positions = work_;
    LatticePosition const origin(0, 0, 0);

    CellContainer::const_iterator i_first = cells_.begin();
    CellContainer::const_iterator const i_end = cells_.end();
    for(; i_first != i_end; ++i_first) {
      positions.resize((*i_first)->GetVertices().size());
      std::fill(positions.begin(), positions.end(), origin);
      velocitiesOnMesh(*i_first, latticeData_, stencil, positions);
      (*i_first)->operator+=(positions);
    }
    // Positions have changed: update Divide and Conquer stuff
    dnc_.update();
  }

template<class KERNEL>
  void CellArmy<KERNEL> :: cell2FluidInteractions() {
    std::vector<LatticeForceVector> &forces = work_;

    CellContainer::const_iterator i_first = cells_.begin();
    CellContainer::const_iterator const i_end = cells_.end();
    for(; i_first != i_end; ++i_first)
      forcesOnGrid<
        typename KERNEL::LatticeType
      >(*i_first, forces, latticeData_, stencil);

    addCell2CellInteractions(dnc_, cell2Cell, stencil, latticeData_);
  }


}}

#endif
