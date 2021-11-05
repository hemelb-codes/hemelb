// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_GRIDANDCELL_H
#define HEMELB_REDBLOOD_GRIDANDCELL_H

#include <vector>
#include <numeric>
#include "units.h"
#include "redblood/Cell.h"
#include "redblood/stencil.h"
#include "redblood/VelocityInterpolation.h"
#include "util/Iterator.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace redblood
  {
// Implementation details
#include "redblood/GridAndCell.impl.h"

    //! Displacement of the cell nodes interpolated from lattice velocities
    template<class KERNEL, class STENCIL>
    void velocitiesOnMesh(std::shared_ptr<CellBase const> cell, geometry::LatticeData const &latDat,
                          std::vector<LatticePosition> &displacements)
    {
      displacements.resize(cell->GetNumberOfNodes());
      std::transform(cell->GetVertices().begin(),
                     cell->GetVertices().end(),
                     displacements.begin(),
                     [&latDat](LatticePosition const &position)
                     {
                       return interpolateVelocity<KERNEL, STENCIL>(latDat, position);
                     });
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
    template<class LATTICE, class STENCIL>
    Dimensionless forcesOnGrid(std::shared_ptr<CellBase const> cell,
                               std::vector<LatticeForceVector> &forces,
                               geometry::LatticeData &latticeData)
    {
      forces.resize(cell->GetNumberOfNodes());
      std::fill(forces.begin(), forces.end(), LatticeForceVector(0, 0, 0));
      auto const energy = cell->Energy(forces);

      typedef details::SpreadForces Spreader;
      details::spreadForce2Grid<Spreader, STENCIL>(cell, Spreader(forces, latticeData));
      return energy;
    }

    //! Computes and Spreads the forces from the cell to the lattice
    //! Adds in the node-wall interaction. It is easier to add here since we
    //! already have a loop over neighboring grid nodes. Assumption is that the
    //! interaction distance is smaller or equal to stencil.
    //! Returns the energy (excluding node-wall interaction)
    template<class LATTICE, class STENCIL>
    Dimensionless forcesOnGrid(std::shared_ptr<CellBase const> cell,
                               geometry::LatticeData &latticeData)
    {
      std::vector<LatticeForceVector> forces(cell->GetNumberOfNodes(), {0.0, 0.0, 0.0});
      return forcesOnGrid<LATTICE, STENCIL>(cell, forces, latticeData);
    }
  }
} // hemelb::redblood

#endif
