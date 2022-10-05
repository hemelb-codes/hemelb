// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "redblood/parallel/GraphBasedCommunication.h"

#include "net/MpiCommunicator.h"
#include "reporting/Timers.h"
#include "util/Iterator.h"

namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
        GlobalCoordsToProcMap ComputeGlobalCoordsToProcMap(net::MpiCommunicator const &comm,
                                                           const geometry::Domain& domain)
        {
          GlobalCoordsToProcMap coordsToProcMap;

          // Populate map with coordinates of locally owned lattice sites first
          std::vector<LatticeVector> locallyOwnedSites;
          locallyOwnedSites.reserve(domain.GetLocalFluidSiteCount());
          auto myRank = comm.Rank();
          for (site_t localSiteId = 0; localSiteId < domain.GetLocalFluidSiteCount(); ++localSiteId)
          {
            auto const& globalSiteCoords = domain.GetSite(localSiteId).GetGlobalSiteCoords();
            locallyOwnedSites.push_back(globalSiteCoords);
            coordsToProcMap[globalSiteCoords] = myRank;
          }

          // Exchange coordinates of locally owned lattice sites with neighbours in comms graph
          auto const& neighbouringProcs = comm.GetNeighbors();
          if (neighbouringProcs.size() > 0)
          {
            auto neighSites = comm.AllNeighGatherV(locallyOwnedSites);
            assert(neighSites.size() == comm.GetNeighborsCount());

            // Finish populating map with knowledge comming from neighbours
            for (auto&& [i, p]: util::enumerate(neighbouringProcs)) {
                for (auto const& globalCoord: neighSites[i]) {
                    // lattice sites are uniquely owned, so no chance of coordinates being repeated across processes
                    assert(coordsToProcMap.count(globalCoord) == 0);
                    coordsToProcMap[globalCoord] = p;
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
                                                                    geometry::Domain &domain,
                                                                    LatticeDistance cellsEffectiveSize)
        {
          std::vector<LatticeVector> serialisedLocalCoords;
          serialisedLocalCoords.reserve(domain.GetDomainEdgeCollisionCount(0));

          for (auto siteIndex = domain.GetMidDomainSiteCount();
               siteIndex < domain.GetMidDomainSiteCount() + domain.GetDomainEdgeCollisionCount(0);
              ++siteIndex)
          {
            serialisedLocalCoords.push_back(domain.GetSite(siteIndex).GetGlobalSiteCoords());
          }

          auto coordsPerProc = comm.AllGatherV(serialisedLocalCoords);

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

          auto const numProcs = comm.Size();
          std::vector<std::vector<int>> vertices(numProcs);
          for (int procA = 0; procA < numProcs; ++procA)
          {
            for (int procB = procA + 1; procB < numProcs; ++procB)
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
                                             geometry::Domain &domain,
                                             std::shared_ptr<TemplateCellContainer> cellTemplates,
                                             hemelb::reporting::Timers &timings)
        {
          timings[hemelb::reporting::Timers::graphComm].Start();
          auto graphComm =
              comm.Graph(ComputeProcessorNeighbourhood(comm,
                                                       domain,
                                                       ComputeCellsEffectiveSize(cellTemplates)));
          timings[hemelb::reporting::Timers::graphComm].Stop();

          return graphComm;
        }

    }
  }
}
