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


namespace hemelb
{
  namespace redblood
  {
    // Implementation details
#include "redblood/GridAndCell.impl.h"

    //! Displacement of the cell nodes interpolated from lattice velocities
    template<class KERNEL> void velocitiesOnMesh(
      CellBase const &cell,
      geometry::LatticeData const &latDat,
      stencil::types stencil,
      std::vector<LatticePosition> &displacements
    )
    {
      displacements.resize(cell.GetNumberOfNodes());
      details::VelocityNodeLoop<KERNEL>(stencil, cell, latDat)
      .loop(details::transform_iterator(displacements.begin()));
    }
    //! Displacement of the cell nodes interpolated from lattice velocities
    template<class KERNEL> void velocitiesOnMesh(
      std::shared_ptr<CellBase const> const &cell,
      geometry::LatticeData const &latDat,
      stencil::types stencil,
      std::vector<LatticePosition> &displacements
    )
    {
      return velocitiesOnMesh<KERNEL>(
               *cell.get(), latDat, stencil, displacements);
    }

    //! \brief Computes and Spreads the forces from the cell to the lattice
    //! \details Adds in the node-wall interaction. It is easier to add here since
    //! already have a loop over neighboring grid nodes. Assumption is that the
    //! interaction distance is smaller or equal to stencil.
    //! \param[in] cell: the cell for which to compute and spread forces
    //! \param[inout] forces: a work array, resized and set to zero prior to use
    //! \param[inout] latticeData: the LB grid
    //! \param[inout] stencil: type of stencil to use when spreading forces
    //! \returns the energy (excluding node-wall interaction)
    template<class LATTICE> Dimensionless forcesOnGrid(
      CellBase const &cell,
      std::vector<LatticeForceVector> &forces,
      geometry::LatticeData &latticeData,
      stencil::types stencil
    )
    {
      forces.resize(cell.GetNumberOfNodes());
      std::fill(forces.begin(), forces.end(), LatticeForceVector(0, 0, 0));
      return cell(forces);
    }

    //! \brief Computes and Spreads the forces from the cell to the lattice
    //! \details Adds in the node-wall interaction. It is easier to add here since
    //! already have a loop over neighboring grid nodes. Assumption is that the
    //! interaction distance is smaller or equal to stencil.
    //! \param[in] cell: the cell for which to compute and spread forces
    //! \param[inout] forces: a work array, resized and set to zero prior to use
    //! \param[inout] latticeData: the LB grid
    //! \param[inout] stencil: type of stencil to use when spreading forces
    //! \returns the energy (excluding node-wall interaction)
    template<class LATTICE> Dimensionless forcesOnGrid(
      Cell const &cell,
      std::vector<LatticeForceVector> &forces,
      geometry::LatticeData &latticeData,
      stencil::types stencil
    )
    {
      Dimensionless const energy = forcesOnGrid<LATTICE>(
                                     *static_cast<CellBase const*>(&cell), forces, latticeData, stencil);

      typedef details::SpreadForcesAndWallForces<LATTICE> Spreader;
      details::spreadForce2Grid(
        cell,
        Spreader(cell, forces, latticeData),
        stencil
      );
      return energy;
    }

    //! \brief Computes and Spreads the forces from the cell to the lattice
    template<class LATTICE> Dimensionless forcesOnGrid(
      std::shared_ptr<Cell const> const &cell,
      std::vector<LatticeForceVector> &forces,
      geometry::LatticeData &latticeData,
      stencil::types stencil
    )
    {
      return forcesOnGrid<LATTICE>(*cell, forces, latticeData, stencil);
    }
    //! \brief Computes and Spreads the forces from the cell to the lattice
    template<class LATTICE> Dimensionless forcesOnGrid(
      std::shared_ptr<CellBase const> const &cell,
      std::vector<LatticeForceVector> &forces,
      geometry::LatticeData &latticeData,
      stencil::types stencil
    )
    {
      return forcesOnGrid<LATTICE>(*cell, forces, latticeData, stencil);
    }

    //! Computes and Spreads the forces from the cell to the lattice
    //! Adds in the node-wall interaction. It is easier to add here since we
    //! already have a loop over neighboring grid nodes. Assumption is that the
    //! interaction distance is smaller or equal to stencil.
    //! Returns the energy (excluding node-wall interaction)
    template<class LATTICE> Dimensionless forcesOnGrid(
      Cell const &cell,
      geometry::LatticeData &latticeData,
      stencil::types stencil
    )
    {
      std::vector<LatticeForceVector> forces(cell.GetNumberOfNodes(), 0);
      return forcesOnGrid<LATTICE>(cell, forces, latticeData, stencil);
    }

  }
} // hemelb::redblood

#endif
