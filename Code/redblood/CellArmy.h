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
#include <iomanip>

#include "redblood/Cell.h"
#include "redblood/CellCell.h"
#include "redblood/WallCellPairIterator.h"
#include "redblood/GridAndCell.h"
#include "redblood/FlowExtension.h"
#include "geometry/LatticeData.h"
#include "redblood/types.h"
#include "redblood/parallel/SpreadForces.h"

namespace hemelb
{
  namespace redblood
  {
    //! \brief Generates a graph communicator describing the data dependencies for interpolation and spreading
    //! @todo Move declaration somewhere more suitable
    //! @todo This is the most conservative and inefficient implementation of the method possible
    net::MpiCommunicator CreateGraphComm(net::MpiCommunicator const &comm)
    {
      if (comm.Size() == 1)
      {
      }
      // setups a graph communicator that in-practice is all-to-all
      // Simpler than setting up something realistic
      std::vector<std::vector<int>> vertices;
      for (size_t i(0); i < comm.Size(); ++i)
      {
        vertices.push_back(std::vector<int>());
        for (size_t j(0); j < comm.Size(); ++j)
        {
          if (j != i)
          {
            vertices[i].push_back(j);
          }
        }
      }
      return comm.Graph(vertices);
    }

    //! \brief Federates the cells together so we can apply ops simultaneously
    //! \tparam TRAITS holds type of kernel and stencil
    template<class TRAITS> class CellArmy
    {
      public:
        typedef typename TRAITS::Lattice Lattice;
        typedef typename TRAITS::Kernel Kernel;
        typedef typename TRAITS::Stencil Stencil;
        //! Type of callback for listening to changes to cells
        typedef std::function<void(const CellContainer &)> CellChangeListener;

        CellArmy(geometry::LatticeData &_latDat, CellContainer const &cells,
                 LatticeDistance boxsize = 10.0, Node2NodeForce const &cell2Cell = { 0e0, 1e0, 2 },
                 Node2NodeForce const &cell2Wall = { 0e0, 1e0, 2 }) :
            latticeData(_latDat), cells(cells), cellDnC(cells, boxsize, cell2Cell.cutoff + 1e-6),
                wallDnC(createWallNodeDnC<Lattice>(_latDat, boxsize, cell2Wall.cutoff + 1e-6)),
                cell2Cell(cell2Cell), cell2Wall(cell2Wall), neighbourDependenciesGraph(CreateGraphComm(net::MpiCommunicator::World()))
        {
        }

        //! Performs fluid to lattice interactions
        void Fluid2CellInteractions();

        //! Performs lattice to fluid interactions
        void Cell2FluidInteractions();

        CellContainer::size_type size()
        {
          return cells.size();
        }

#   ifdef HEMELB_DOING_UNITTESTS
        //! Updates divide and conquer
        void updateDNC()
        {
          cellDnC.update();
        }
        CellContainer const & GetCells() const
        {
          return cells;
        }
        DivideConquerCells const & GetDNC() const
        {
          return cellDnC;
        }
#   endif

        //! Sets up call for cell insertion
        //! Called everytime CallCellInsertion is called
        void SetCellInsertion(std::function<void(CellInserter const&)> const & f)
        {
          cellInsertionCallBack = f;
        }

        //! Calls cell insertion
        void CallCellInsertion()
        {
          if (cellInsertionCallBack)
          {
            auto callback = [this](CellContainer::value_type cell)
            {
              this->AddCell(cell);
            };
            cellInsertionCallBack(callback);
          }
        }

        //! Adds a cell change listener to be notified when cell positions change
        void AddCellChangeListener(CellChangeListener const & ccl)
        {
          cellChangeListeners.push_back(ccl);
        }

        //! Invokes the callback function to output cell positions
        void NotifyCellChangeListeners()
        {
          for (CellChangeListener ccl : cellChangeListeners)
          {
            ccl(cells);
          }
        }

        //! Sets outlets within which cells disappear
        void SetOutlets(std::vector<FlowExtension> const & olets)
        {
          outlets = olets;
        }

        //! Remove cells if they have reached outlets
        void CellRemoval();

        //! Adds input cell to simulation
        void AddCell(CellContainer::value_type cell)
        {
          auto const barycenter = cell->GetBarycenter();
          log::Logger::Log<log::Debug, log::OnePerCore>("Adding cell at (%f, %f, %f)",
                                                        barycenter.x,
                                                        barycenter.y,
                                                        barycenter.z);
          cellDnC.insert(cell);
          cells.insert(cell);
        }

        //! \brief Sets cell to cell interaction forces
        //! \details Forwards arguments to Node2NodeForce constructor.
        template<class ... ARGS> void SetCell2Cell(ARGS && ... args)
        {
          cell2Cell = Node2NodeForce(std::forward<ARGS>(args)...);
          cellDnC.SetBoxSizeAndHalo(cellDnC.GetBoxSize(), cell2Cell.cutoff + 1e-6);
        }
        //! \brief Sets cell to cell interaction forces
        //! \details Forwards arguments to Node2NodeForce constructor.
        template<class ... ARGS> void SetCell2Wall(ARGS && ... args)
        {
          cell2Wall = Node2NodeForce(std::forward<ARGS>(args)...);
          wallDnC = createWallNodeDnC<Lattice>(wallDnC,
                                               wallDnC.GetBoxSize(),
                                               cell2Wall.cutoff + 1e-6);
        }

      protected:
        //! All lattice information and then some
        geometry::LatticeData &latticeData;
        //! Contains all cells
        CellContainer cells;
        //! Divide and conquer object
        DivideConquerCells cellDnC;
        //! Divide and conquer object
        DivideConquer<WallNode> wallDnC;
        //! A work array with forces/positions
        std::vector<LatticePosition> work;
        //! This function is called every lb turn
        //! It should insert cells using the call back passed to it.
        std::function<void(CellInserter const&)> cellInsertionCallBack;
        //! Observers to be notified when cell positions change
        std::vector<CellChangeListener> cellChangeListeners;
        //! Remove cells if they reach these outlets
        std::vector<FlowExtension> outlets;
        //! Interaction terms between cells
        Node2NodeForce cell2Cell;
        //! Interaction terms between cells
        Node2NodeForce cell2Wall;
        //! Communicator defining the data depencies between processors for spreading/interpolation
        net::MpiCommunicator neighbourDependenciesGraph;

    };

    template<class TRAITS>
    void CellArmy<TRAITS>::Fluid2CellInteractions()
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Fluid -> cell interations");

      auto const distributions = hemelb::redblood::parallel::nodeDistributions(latticeData, cells);

      ExchangeCells xc(neighbourDependenciesGraph);
      auto ownership = [this](CellContainer::value_type cell) {
        auto const id = latticeData.GetProcIdFromGlobalCoords(cell->GetBarycenter());
        if (id != BIG_NUMBER2)
        {
          throw Exception("Cell nobody owns");
        }
        return id;
      };
      xc.PostCellMessageLength(distributions, cells, ownership);
      xc.PostCells(distributions, cells, ownership);
      auto const distCells = xc.ReceiveCells(rbcMeshes);
      xc.Update(cells, result);
      auto const &lentCells = std::get<2>(distCells);

      // Actually perform velocity integration
      hemelb::redblood::parallel::IntegrateVelocities integrator(neighbourDependenciesGraph);
      integrator.PostMessageLength(lentCells);
      integrator.ComputeLocalVelocitiesAndUpdatePositions<Traits>(latticeData, cells);
      integrator.PostVelocities<Traits>(latticeData, lentCells);
      integrator.UpdatePositionsNonLocal(distributions, cells);

      //! @todo Any changes required for these lines when running in parallel?
      // Positions have changed: update Divide and Conquer stuff
      // update should be a function taking the distCells tuple (newly owned, newly disowned, lent)
      // and update the cell container within cellDnC
      cellDnC.update();
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::Cell2FluidInteractions()
    {
      using namespace hemelb::redblood::parallel;

      log::Logger::Log<log::Debug, log::OnePerCore>("Cell -> fluid interations");
      latticeData.ResetForces();

      auto const distributions = nodeDistributions(latticeData, cells);

      SpreadForces mpi_spreader(neighbourDependenciesGraph);
      mpi_spreader.PostMessageLength(distributions, cells);
      mpi_spreader.ComputeForces(cells);
      mpi_spreader.PostForcesAndNodes(distributions, cells);
      mpi_spreader.SpreadLocalForces<TRAITS>(latticeData, cells);
      mpi_spreader.SpreadNonLocalForces<TRAITS>(latticeData);

      //! @todo Any changes required for these lines when running in parallel?
      addCell2CellInteractions<Stencil>(cellDnC, cell2Cell, latticeData);
      addCell2WallInteractions<Stencil>(cellDnC, wallDnC, cell2Wall, latticeData);
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::CellRemoval()
    {
      auto i_first = cells.cbegin();
      auto const i_end = cells.cend();
      while (i_first != i_end)
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
        if (std::find_if(outlets.begin(), outlets.end(), checkCell) != outlets.end())
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("Removing cell at (%f, %f, %f)",
                                                        barycenter.x,
                                                        barycenter.y,
                                                        barycenter.z);
          cellDnC.remove(*i_current);
          cells.erase(i_current);
        }
      }
    }
  }
}

#endif
