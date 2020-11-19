// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "redblood/parallel/GraphBasedCommunication.h"
#include "util/Iterator.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
        GlobalCoordsToProcMap ComputeGlobalCoordsToProcMap(net::MpiCommunicator const &comm,
                                                           const geometry::LatticeData &latDat)
        {
          GlobalCoordsToProcMap coordsToProcMap;

          // Populate map with coordinates of locally owned lattice sites first
          std::vector<LatticeVector> locallyOwnedSites;
          locallyOwnedSites.reserve(latDat.GetLocalFluidSiteCount());
          auto myRank = comm.Rank();
          for (site_t localSiteId = 0; localSiteId < latDat.GetLocalFluidSiteCount(); ++localSiteId)
          {
            auto const& globalSiteCoords = latDat.GetSite(localSiteId).GetGlobalSiteCoords();
            locallyOwnedSites.push_back(globalSiteCoords);
            coordsToProcMap[globalSiteCoords] = myRank;
          }

          // Exchange coordinates of locally owned lattice sites with neighbours in comms graph
          auto const& neighbouringProcs = comm.GetNeighbors();
          if (neighbouringProcs.size() > 0)
          {
            std::vector<std::vector<LatticeVector>> neighSites = comm.AllNeighGatherV(locallyOwnedSites);
            assert(neighSites.size() == comm.GetNeighborsCount());

            // Finish populating map with knowledge comming from neighbours
            for (auto const& neighbour : hemelb::util::enumerate(neighbouringProcs))
            {
              for (auto const& globalCoord : neighSites[neighbour.index])
              {
                // lattice sites are uniquely owned, so no chance of coordinates being repeated across processes
                assert(coordsToProcMap.count(globalCoord) == 0);
                coordsToProcMap[globalCoord] = neighbour.value;
              }
            }
          }

          return coordsToProcMap;
        }

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
          std::vector<LatticeVector> serialisedLocalCoords;
          serialisedLocalCoords.reserve(latDat.GetDomainEdgeCollisionCount(0));

          for (auto siteIndex = latDat.GetMidDomainSiteCount();
              siteIndex < latDat.GetMidDomainSiteCount() + latDat.GetDomainEdgeCollisionCount(0);
              ++siteIndex)
          {
            serialisedLocalCoords.push_back(latDat.GetSite(siteIndex).GetGlobalSiteCoords());
          }

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

          std::vector<LatticeVector> allSerialisedCoords(totalSize);
          HEMELB_MPI_CALL(MPI_Allgatherv,
                          ( net::MpiConstCast(&serialisedLocalCoords[0]), serialisedLocalCoords.size(), net::MpiDataType<LatticeVector>(), &allSerialisedCoords[0], net::MpiConstCast(&allSerialisedCoordSizes[0]), net::MpiConstCast(&allSerialisedCoordDisplacements[0]), net::MpiDataType<LatticeVector>(), comm ));

          std::vector<std::vector<LatticeVector>> coordsPerProc(numProcs);
          for (decltype(numProcs) procIndex = 0; procIndex < numProcs; ++procIndex)
          {
            for (auto indexAllCoords = allSerialisedCoordDisplacements[procIndex];
                 indexAllCoords < allSerialisedCoordDisplacements[procIndex + 1]; ++indexAllCoords)
            {
              coordsPerProc[procIndex].push_back(allSerialisedCoords[indexAllCoords]);
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
            for (int procB(procA+1); procB < numProcs; ++procB)
            {
              if (areProcsNeighbours(procA, procB))
              {
                vertices[procA].push_back(procB);
                vertices[procB].push_back(procA);
              }
            }
          }

          return vertices;
        }

        LatticeDistance ComputeCellsEffectiveSize(std::shared_ptr<TemplateCellContainer> cellTemplates)
        {
          double maxCellRadius = std::numeric_limits<LatticeDistance>::min();

          for (auto cellTemplate : *cellTemplates)
          {
            maxCellRadius = std::max(maxCellRadius, cellTemplate.second->GetScale());
          }

          return MAXIMUM_SIZE_TO_RADIUS_RATIO * maxCellRadius;
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

    }
  }
}
