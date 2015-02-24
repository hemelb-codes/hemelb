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
#include "redblood/GridAndCell.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace redblood
  {

    //! \brief Federates the cells together so we can apply ops simultaneously
    template<class KERNEL> class CellArmy
    {
      public:
        //! Interaction terms between cells
        Node2NodeForce cell2Cell;
        //! Stencil
        stencil::types stencil = stencil::types::FOUR_POINT;

        CellArmy(geometry::LatticeData &_latDat,
            CellContainer const &cells,
            PhysicalDistance boxsize=10.0, PhysicalDistance halo=2.0)
          : latticeData(_latDat), cells(cells),
          dnc(cells, boxsize, halo)
        {
        }

        //! Performs fluid to lattice interactions
        void Fluid2CellInteractions();

        //! Performs lattice to fluid interactions
        void Cell2FluidInteractions();

    #   ifdef HEMELB_DOING_UNITTESTS
          //! Updates divide and conquer
          void updateDNC()
          {
            dnc.update();
          }
    #   endif

      protected:
        //! All lattice information and then some
        geometry::LatticeData &latticeData;
        //! Contains all cells
        CellContainer cells;
        //! Divide and conquer object
        DivideConquerCells dnc;
        //! A work array with forces/positions
        std::vector<LatticePosition> work;
    };

    template<class KERNEL>
      void CellArmy<KERNEL> :: Fluid2CellInteractions()
      {
        std::vector<LatticePosition> & positions = work;
        LatticePosition const origin(0, 0, 0);

        CellContainer::const_iterator i_first = cells.begin();
        CellContainer::const_iterator const i_end = cells.end();
        for(; i_first != i_end; ++i_first)
        {
          positions.resize((*i_first)->GetVertices().size());
          std::fill(positions.begin(), positions.end(), origin);
          velocitiesOnMesh<KERNEL>(*i_first, latticeData, stencil, positions);
          (*i_first)->operator+=(positions);
        }
        // Positions have changed: update Divide and Conquer stuff
        dnc.update();
      }

    template<class KERNEL>
      void CellArmy<KERNEL> :: Cell2FluidInteractions()
      {
        std::vector<LatticeForceVector> &forces = work;

        CellContainer::const_iterator i_first = cells.begin();
        CellContainer::const_iterator const i_end = cells.end();
        for(; i_first != i_end; ++i_first)
        {
          forcesOnGrid<typename KERNEL::LatticeType>(*i_first, forces, latticeData, stencil);
        }

        addCell2CellInteractions(dnc, cell2Cell, stencil, latticeData);
      }


}}

#endif
