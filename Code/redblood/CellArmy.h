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
#include "redblood/FlowExtension.h"
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

        //! Outputs cell positions
        void CellOutput(std::ostream &) const;

    #   ifdef HEMELB_DOING_UNITTESTS
          //! Updates divide and conquer
          void updateDNC()
          {
            dnc.update();
          }
          CellContainer const & GetCells() const
          {
            return cells;
          }
          DivideConquerCells const & GetDNC() const
          {
            return dnc;
          }
    #   endif

        //! Sets up call for cell insertion
        //! Called everytime CallCellInsertion is called
        void SetCellInsertion(std::function<void(std::function<void(CellContainer::value_type)>)> f)
        {
          cellInsertionCallBack = f;
        }

        //! Calls cell insertion
        void CallCellInsertion()
        {
          if(cellInsertionCallBack)
            cellInsertionCallBack([this](CellContainer::value_type cell) { this->AddCell(cell); });
        }

        //! Sets outlets within which cells disapear
        void SetOutlets(std::vector<FlowExtension> const & olets)
        {
          outlets = olets;
        }

        //! Remove cells if they have reached outlets
        void CellRemoval();

      protected:
        //! Adds input cell to simulation
        void AddCell(CellContainer::value_type cell)
        {
          dnc.insert(cell);
          cells.insert(cell);
        }

        //! All lattice information and then some
        geometry::LatticeData &latticeData;
        //! Contains all cells
        CellContainer cells;
        //! Divide and conquer object
        DivideConquerCells dnc;
        //! A work array with forces/positions
        std::vector<LatticePosition> work;
        //! This function is called every lb turn
        //! It should insert cells using the call back passed to it.
        std::function<void(std::function<void(CellContainer::value_type)>)>
          cellInsertionCallBack;

        //! Remove cells if they reach these outlets
        std::vector<FlowExtension> outlets;
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

    template<class KERNEL>
      void CellArmy<KERNEL> :: CellOutput(std::ostream & out) const
      {
        CellContainer::const_iterator end = cells.cend();
        for (CellContainer::const_iterator i = cells.cbegin(); i != end; ++i)
        {
          MeshData::Vertices::value_type centre = (*i)->GetBarycenter();
          out << centre << std::endl;
        }
      }

    template<class KERNEL>
      void CellArmy<KERNEL> :: CellRemoval()
      {
        auto i_first = cells.cbegin();
        auto const i_end = cells.cend();
        while(i_first != i_end)
        {
          auto const barycenter = (*i_first)->GetBarycenter();
          auto checkCell = [&barycenter](FlowExtension const &flow)
          {
            return contains(flow, barycenter);
          };
          // save current iterator and increment before potential removal.
          // removing the cell from the set should invalidate only the relevant iterator.
          auto const i_current = i_first;
          ++i_first;
          if(std::find_if(outlets.begin(), outlets.end(), checkCell) != outlets.end())
          {
            dnc.remove(*i_current);
            cells.erase(i_current);
          }
        }
      }
}}

#endif
