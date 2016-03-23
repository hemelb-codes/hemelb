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
#include "Exception.h"
#include "redblood/parallel/IntegrateVelocities.h"

namespace hemelb
{
  namespace redblood
  {
    //! \brief All processes are considered neighbours with each other. This is the most conservative and inefficient implementation of the method possible.
    std::vector<std::vector<int>> ComputeProcessorNeighbourhood(net::MpiCommunicator const &comm)
    {
      // setups a graph communicator that in-practice is all-to-all
      std::vector<std::vector<int>> vertices;
      for (int i(0); i < comm.Size(); ++i)
      {
        vertices.push_back(std::vector<int>());
        for (int j(0); j < comm.Size(); ++j)
        {
          if (j != i)
          {
            vertices[i].push_back(j);
          }
        }
      }

      return vertices;
    }

    //! \brief Compute neighbourhood based on checking the minimum distance between every pair of subdomains and declaring them neighbours if this is shorter than the RBCs effective size.
    std::vector<std::vector<int>> ComputeProcessorNeighbourhood(net::MpiCommunicator const &comm,
                                                                geometry::LatticeData &latDat,
                                                                LatticeDistance cellsEffectiveSize)
    {
      std::vector<LatticeCoordinate> serialisedLocalCoords(3
          * latDat.GetDomainEdgeCollisionCount(0));
      std::vector<LatticeCoordinate>::size_type serialisedLocalCoordsIndex = 0;

      for (auto siteIndex = latDat.GetMidDomainSiteCount();
          siteIndex < latDat.GetMidDomainSiteCount() + latDat.GetDomainEdgeCollisionCount(0);
          ++siteIndex)
      {
        auto borderSiteCoords = latDat.GetSite(siteIndex).GetGlobalSiteCoords();
        serialisedLocalCoords[serialisedLocalCoordsIndex++] = borderSiteCoords[0];
        serialisedLocalCoords[serialisedLocalCoordsIndex++] = borderSiteCoords[1];
        serialisedLocalCoords[serialisedLocalCoordsIndex++] = borderSiteCoords[2];
      }
      assert(serialisedLocalCoordsIndex == serialisedLocalCoords.size());

      /// @\todo refactor into a method net::MpiCommunicator::AllGatherv
      int numProcs = comm.Size();
      std::vector<int> allSerialisedCoordSizes = comm.AllGather((int) serialisedLocalCoords.size());
      std::vector<int> allSerialisedCoordDisplacements(numProcs + 1);

      site_t totalSize = std::accumulate(allSerialisedCoordSizes.begin(),
                                         allSerialisedCoordSizes.end(),
                                         0);

      allSerialisedCoordDisplacements[0] = 0;
      for (int j = 0; j < numProcs; ++j)
      {
        allSerialisedCoordDisplacements[j + 1] = allSerialisedCoordDisplacements[j]
            + allSerialisedCoordSizes[j];
      }

      std::vector<LatticeCoordinate> allSerialisedCoords(totalSize);
      HEMELB_MPI_CALL(MPI_Allgatherv,
                      ( net::MpiConstCast(&serialisedLocalCoords[0]), serialisedLocalCoords.size(), net::MpiDataType<LatticeCoordinate>(), &allSerialisedCoords[0], net::MpiConstCast(&allSerialisedCoordSizes[0]), net::MpiConstCast(&allSerialisedCoordDisplacements[0]), net::MpiDataType<LatticeCoordinate>(), comm ));

      std::vector<std::vector<LatticeVector>> coordsPerProc(numProcs);
      for (unsigned procIndex = 0; procIndex < numProcs; ++procIndex)
      {
        for (unsigned indexAllCoords = allSerialisedCoordDisplacements[procIndex];
            indexAllCoords < allSerialisedCoordDisplacements[procIndex + 1]; indexAllCoords += 3)
        {
          coordsPerProc[procIndex].push_back( { allSerialisedCoords[indexAllCoords],
                                                allSerialisedCoords[indexAllCoords + 1],
                                                allSerialisedCoords[indexAllCoords + 2] });
        }
      }
      /// end of refactoring

      auto cellsEffectiveSizeSq = cellsEffectiveSize * cellsEffectiveSize;
      auto areProcsNeighbours =
          [cellsEffectiveSizeSq, &coordsPerProc] (unsigned procA, unsigned procB)
          {
            if (procA == procB)
            {
              return false;
            }

            auto distanceSqBetweenSubdomainEdges = std::numeric_limits<LatticeDistance>::max();
            for(auto siteProcA : coordsPerProc[procA])
            {
              for(auto siteProcB : coordsPerProc[procB])
              {
                distanceSqBetweenSubdomainEdges = std::min(distanceSqBetweenSubdomainEdges, (LatticeDistance)(siteProcA-siteProcB).GetMagnitudeSquared());
              }
            }

            return distanceSqBetweenSubdomainEdges < cellsEffectiveSizeSq;
          };

      std::vector<std::vector<int>> vertices(numProcs);
      for (int procA(0); procA < numProcs; ++procA)
      {
        for (int procB(0); procB < numProcs; ++procB)
        {
          if (areProcsNeighbours(procA, procB))
          {
            vertices[procA].push_back(procB);
          }
        }
      }

      return vertices;
    }

    // Make effective size 1.5 times the diameter
    static const LatticeDistance EFFECTIVE_SIZE_TO_RADIUS_RATIO = 3.0;

    LatticeDistance ComputeCellsEffectiveSize(std::shared_ptr<TemplateCellContainer> cellTemplates)
    {
      double maxCellRadius = std::numeric_limits<LatticeDistance>::min();

      for (auto cellTemplate : *cellTemplates)
      {
        maxCellRadius = std::max(maxCellRadius, cellTemplate.second->GetScale());
      }

      return EFFECTIVE_SIZE_TO_RADIUS_RATIO * maxCellRadius;
    }

    //! \brief Generates a graph communicator describing the data dependencies for interpolation and spreading
    net::MpiCommunicator CreateGraphComm(net::MpiCommunicator const &comm,
                                         geometry::LatticeData &latDat,
                                         std::shared_ptr<TemplateCellContainer> cellTemplates,
                                         hemelb::reporting::Timers &timings)
    {
      timings[hemelb::reporting::Timers::graphComm].Start();
      auto graphComm =
          comm.Graph(ComputeProcessorNeighbourhood(comm,
                                                   latDat,
                                                   ComputeCellsEffectiveSize(cellTemplates)));
      timings[hemelb::reporting::Timers::graphComm].Stop();

      return graphComm;
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

        CellArmy(geometry::LatticeData &latDat, CellContainer const &cells,
                 std::shared_ptr<TemplateCellContainer> cellTemplates,
                 hemelb::reporting::Timers &timings, LatticeDistance boxsize = 10.0,
                 Node2NodeForce const &cell2Cell = { 0e0, 1e0, 2 },
                 Node2NodeForce const &cell2Wall = { 0e0, 1e0, 2 },
                 net::MpiCommunicator const &worldCommunicator = net::MpiCommunicator::World()) :
            latticeData(latDat), cells(cells), cellDnC(cells, boxsize, cell2Cell.cutoff + 1e-6),
                wallDnC(createWallNodeDnC<Lattice>(latDat, boxsize, cell2Wall.cutoff + 1e-6)),
                cell2Cell(cell2Cell), cell2Wall(cell2Wall), worldCommunicator(worldCommunicator),
                cellTemplates(cellTemplates), timings(timings),
                neighbourDependenciesGraph(CreateGraphComm(worldCommunicator,
                                                           latDat,
                                                           cellTemplates,
                                                           timings)),
                exchangeCells(neighbourDependenciesGraph, worldCommunicator),
                velocityIntegrator(neighbourDependenciesGraph),
                forceSpreader(neighbourDependenciesGraph),
                nodeDistributions(parallel::nodeDistributions(latticeData, cells))
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
        parallel::CellParallelization::LentCells const & GetLentCells() const
        {
          return lentCells;
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

        std::function<void(CellInserter const&)> GetCellInsertion() const
        {
          return cellInsertionCallBack;
        }

        //! Calls cell insertion
        void CallCellInsertion()
        {
          timings[hemelb::reporting::Timers::cellInsertion].Start();
          if (cellInsertionCallBack)
          {
            auto callback = [this](CellContainer::value_type cell)
            {
              this->AddCell(cell);
            };
            cellInsertionCallBack(callback);
          }
          timings[hemelb::reporting::Timers::cellInsertion].Stop();
        }

        //! Adds a cell change listener to be notified when cell positions change
        void AddCellChangeListener(CellChangeListener const & ccl)
        {
          cellChangeListeners.push_back(ccl);
        }

        //! Invokes the callback function to output cell positions
        void NotifyCellChangeListeners()
        {
          timings[hemelb::reporting::Timers::cellListeners].Start();
          for (CellChangeListener ccl : cellChangeListeners)
          {
            ccl(cells);
          }
          timings[hemelb::reporting::Timers::cellListeners].Stop();
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

          //! @todo: #623 AddCell should only be called if the subdomain contains the relevant RBC inlet
          auto const id = latticeData.GetProcIdFromGlobalCoords(barycenter);
          if (id == latticeData.GetCommunicator().Rank())
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Adding cell at (%f, %f, %f)",
                                                          barycenter.x,
                                                          barycenter.y,
                                                          barycenter.z);
            cellDnC.insert(cell);
            cells.insert(cell);

            nodeDistributions.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(cell->GetTag()),
                                      std::forward_as_tuple(parallel::details::AssessMPIFunction<
                                                                Stencil>(latticeData),
                                                            cell));
          }
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
        //! Cells lent to this process
        parallel::CellParallelization::LentCells lentCells;
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
        //! Communicator with all the processes participating in the simulation
        net::MpiCommunicator const &worldCommunicator;
        //! Container with the templates used to create the RBC meshes
        std::shared_ptr<TemplateCellContainer> cellTemplates;
        //! Timers object used to time different code sections
        hemelb::reporting::Timers &timings;
        //! Communicator defining the data dependencies between processors for spreading/interpolation
        net::MpiCommunicator neighbourDependenciesGraph;
        //! Exchange cells object
        parallel::ExchangeCells exchangeCells;
        //! Velocity integrator object
        parallel::IntegrateVelocities velocityIntegrator;
        //! Force spreader object
        parallel::SpreadForces forceSpreader;
        //! Object describing how the cells affect different subdomains
        parallel::CellParallelization::NodeDistributions nodeDistributions;

    };

    template<class TRAITS>
    void CellArmy<TRAITS>::Fluid2CellInteractions()
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Fluid -> cell interations");

      timings[hemelb::reporting::Timers::exchangeCells].Start();
      auto ownership = [this](CellContainer::value_type cell)
      {
        auto const id = latticeData.GetProcIdFromGlobalCoords(cell->GetBarycenter());
        if (id == BIG_NUMBER2)
        {
          throw Exception() << "Cell nobody owns";
        }
        return id;
      };
      exchangeCells.PostCellMessageLength(nodeDistributions, cells, ownership);
      exchangeCells.PostCells(nodeDistributions, cells, ownership);
      auto const distCells = exchangeCells.ReceiveCells(cellTemplates);
      exchangeCells.Update(cells, distCells);
      exchangeCells.Update(nodeDistributions,
                           distCells,
                           parallel::details::AssessMPIFunction<Stencil>(latticeData));
      timings[hemelb::reporting::Timers::exchangeCells].Stop();

      // Actually perform velocity integration
      timings[hemelb::reporting::Timers::computeAndPostVelocities].Start();
      velocityIntegrator.PostMessageLength(std::get<2>(distCells));
      velocityIntegrator.ComputeLocalVelocitiesAndUpdatePositions<TRAITS>(latticeData, cells);
      velocityIntegrator.PostVelocities<TRAITS>(latticeData, std::get<2>(distCells));
      timings[hemelb::reporting::Timers::computeAndPostVelocities].Stop();

      timings[hemelb::reporting::Timers::receiveVelocitiesAndUpdate].Start();
      velocityIntegrator.UpdatePositionsNonLocal(nodeDistributions, cells);
      timings[hemelb::reporting::Timers::receiveVelocitiesAndUpdate].Stop();

      // Positions have changed: update node distributions
      timings[hemelb::reporting::Timers::computeNodeDistributions].Start();
      for (auto cell : cells)
      {
        nodeDistributions.at(cell->GetTag()).template Reindex<Stencil>(latticeData, cell);
      }
      timings[hemelb::reporting::Timers::computeNodeDistributions].Stop();

      // Positions have changed: update Divide and Conquer stuff
      log::Logger::Log<log::Debug, log::OnePerCore>("Number of lent cells: %i",
                                                    std::get<2>(distCells).size());
      timings[hemelb::reporting::Timers::updateDNC].Start();
      cellDnC.update(distCells);
      timings[hemelb::reporting::Timers::updateDNC].Stop();
      lentCells = std::move(std::get<2>(distCells));
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::Cell2FluidInteractions()
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Cell -> fluid interations");
      latticeData.ResetForces();

      timings[hemelb::reporting::Timers::computeAndPostForces].Start();
      forceSpreader.PostMessageLength(nodeDistributions, cells);
      forceSpreader.ComputeForces(cells);
      forceSpreader.PostForcesAndNodes(nodeDistributions, cells);
      forceSpreader.SpreadLocalForces<TRAITS>(latticeData, cells);
      timings[hemelb::reporting::Timers::computeAndPostForces].Stop();

      timings[hemelb::reporting::Timers::receiveForcesAndUpdate].Start();
      forceSpreader.SpreadNonLocalForces<TRAITS>(latticeData);
      timings[hemelb::reporting::Timers::receiveForcesAndUpdate].Stop();

      //! @todo Any changes required for these lines when running in parallel?
      timings[hemelb::reporting::Timers::updateCellAndWallInteractions].Start();
      addCell2CellInteractions<Stencil>(cellDnC, cell2Cell, latticeData);
      addCell2WallInteractions<Stencil>(cellDnC, wallDnC, cell2Wall, latticeData);
      timings[hemelb::reporting::Timers::updateCellAndWallInteractions].Stop();
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::CellRemoval()
    {
      timings[hemelb::reporting::Timers::cellRemoval].Start();
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
      timings[hemelb::reporting::Timers::cellRemoval].Stop();
    }
  }
}

#endif
