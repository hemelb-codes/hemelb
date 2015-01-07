//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_GRID_AND_CELL_H
#define HEMELB_REDBLOOD_GRID_AND_CELL_H

#include <vector>
#include "units.h"
#include "redblood/Cell.h"
#include "redblood/stencil.h"
#include "redblood/VelocityInterpolation.h"
#include "geometry/LatticeData.h"


namespace hemelb { namespace redblood {
// Implementation details
#include "redblood/GridAndCell.impl.h"

//! Displacement of the cell nodes interpolated from lattice velocities
template<class T_KERNEL> void velocitiesOnMesh(
    Cell const &_cell,
    geometry::LatticeData const &_latDat,
    stencil::types _stencil,
    std::vector<LatticePosition> &_displacements
) {
  _displacements.resize(_cell.GetNumberOfNodes());
  details::VelocityNodeLoop<T_KERNEL>(_stencil, _cell, _latDat)
    .loop(details::transform_iterator(_displacements.begin()));
}

//! Computes and Spreads the forces from the cell to the lattice
//! Returns the energy
Dimensionless forcesOnGrid(
    Cell const &_cell,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
);

//! \brief Computes and Spreads the forces from the cell to the lattice
//! \details Adds in the node-wall interaction. It is easier to add here since
//! already have a loop over neighboring grid nodes. Assumption is that the
//! interaction distance is smaller or equal to stencil.
//! \param[in] _cell: the cell for which to compute and spread forces
//! \param[inout] _forces: a work array, resized and set to zero prior to use
//! \param[inout] _latticeData: the LB grid
//! \param[inout] _stencil: type of stencil to use when spreading forces
//! \returns the energy (excluding node-wall interaction)
template<class LATTICE> Dimensionless forcesOnGridWithWallInteraction(
    Cell const &_cell,
    std::vector<LatticeForceVector> &_forces,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
) {
  typedef details::SpreadForcesAndWallForces<LATTICE> Spreader;
  _forces.resize(_cell.GetNumberOfNodes());
  std::fill(_forces.begin(), _forces.end(), LatticeForceVector(0, 0, 0));
  Dimensionless const energy = _cell(_forces);

  details::spreadForce2Grid(
      _cell,
      Spreader(_cell, _forces, _latticeData),
      _stencil
  );
  return energy;
}

//! Computes and Spreads the forces from the cell to the lattice
//! Adds in the node-wall interaction. It is easier to add here since we
//! already have a loop over neighboring grid nodes. Assumption is that the
//! interaction distance is smaller or equal to stencil.
//! Returns the energy (excluding node-wall interaction)
template<class LATTICE> Dimensionless forcesOnGridWithWallInteraction(
    Cell const &_cell,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
) {
  std::vector<LatticeForceVector> forces(_cell.GetNumberOfNodes(), 0);
  return forcesOnGridWithWallInteraction<LATTICE>(
      _cell, forces, _latticeData, _stencil
  );
}

}} // hemelb::redblood

#endif
