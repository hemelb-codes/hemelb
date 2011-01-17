#include "topology/TopologyManager.h"

#include <math.h>

namespace hemelb
{
  namespace topology
  {
    void TopologyManager::DecomposeDomain(int iTotalFluidSites,
                                          NetworkTopology* bNetTop,
                                          const lb::GlobalLatticeData & bGlobLatDat)
    {
      // Allocations.  fluid sites will store actual number of fluid
      // sites per proc.  Site location will store up to 10000 of some
      // sort of coordinate.

      // Initialise the count of fluid sites on each processor to 0.
      bNetTop->FluidSitesOnEachProcessor
          = new int[bNetTop->GetProcessorCount()];

      for (int n = 0; n < bNetTop->GetProcessorCount(); n++)
      {
        bNetTop->FluidSitesOnEachProcessor[n] = 0;
      }

      // Count of fluid sites not yet visited.
      int lUnvisitedFluidSiteCount = iTotalFluidSites;

      // If one machine or one machine per proc.
      if (bNetTop->MachineCount == 1 || bNetTop->MachineCount
          == bNetTop->GetProcessorCount())
      {
        // Fluid sites per rank.
        int fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
            / (double) bNetTop->GetProcessorCount());

        //Rank we're looking at.
        int proc_count = 0;

        // If we're steering with more than one processor, save one processor for doing that.
#ifndef NO_STEER
        if (bNetTop->GetProcessorCount() != 1)
        {
          fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
              / (double) (bNetTop->GetProcessorCount() - 1));
          proc_count = 1;
        }
#endif

        // In the simple case, simply divide fluid sites up between processors.
        AssignFluidSitesToProcessors(proc_count, fluid_sites_per_unit,
                                     lUnvisitedFluidSiteCount, -1, false,
                                     bGlobLatDat, bNetTop);
      }
      else
      {
        // Rank we are looking at.
        int proc_count = bNetTop->GetProcessorCount();
        double weight = (double) (bNetTop->ProcCountOnEachMachine[0]
            * bNetTop->GetProcessorCount())
            / (double) (bNetTop->GetProcessorCount() - 1);
        // Fluid sites per rank.
        int fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
            * weight / bNetTop->MachineCount);

        // First, divide the sites up between machines.
        AssignFluidSitesToProcessors(proc_count, fluid_sites_per_unit,
                                     lUnvisitedFluidSiteCount, -1, true,
                                     bGlobLatDat, bNetTop);

        fluid_sites_per_unit = (int) ceil((double) lUnvisitedFluidSiteCount
            / (double) (bNetTop->GetProcessorCount() - 1));
        proc_count = 1;

        // For each machine, divide up the sites it has between its cores.
        for (int lMachineNumber = 0; lMachineNumber < bNetTop->MachineCount; lMachineNumber++)
        {
          AssignFluidSitesToProcessors(proc_count, fluid_sites_per_unit,
                                       lUnvisitedFluidSiteCount,
                                       bNetTop->GetProcessorCount()
                                           + lMachineNumber, false,
                                       bGlobLatDat, bNetTop);
        }
      }

    }

    void TopologyManager::AssignFluidSitesToProcessors(int & proc_count,
                                                       int & iSitesPerProc,
                                                       int & bUnassignedSites,
                                                       const int iMarker,
                                                       const bool iIsMachineLevel,
                                                       const lb::GlobalLatticeData &iGlobLatDat,
                                                       NetworkTopology * bNetTopology)
    {
      std::vector<SiteLocation*> *lSiteLocationA = new std::vector<
          SiteLocation*>;
      std::vector<SiteLocation*> *lSiteLocationB = new std::vector<
          SiteLocation*>;

      int lBlockNumber = -1;

      // Domain Decomposition.  Pick a site. Set it to the rank we are
      // looking at. Find its neighbours and put those on the same
      // rank, then find the next-nearest neighbours, etc. until we
      // have a completely joined region, or there are enough fluid
      // sites on the rank.  In the former case, start again at
      // another site. In the latter case, move on to the next rank.
      // Do this until all sites are assigned to a rank. There is a
      // high chance of of all sites on a rank being joined.

      int lSitesOnCurrentProc = 0;

      // Iterate over all blocks.
      for (int lBlockCoordI = 0; lBlockCoordI < iGlobLatDat.GetXBlockCount(); lBlockCoordI++)
      {
        for (int lBlockCoordJ = 0; lBlockCoordJ < iGlobLatDat.GetYBlockCount(); lBlockCoordJ++)
        {
          for (int lBlockCoordK = 0; lBlockCoordK
              < iGlobLatDat.GetZBlockCount(); lBlockCoordK++)
          {
            // Block number is the number of the block we're currently on.
            lBlockNumber++;

            // Point to a block of ProcessorRankForEachBlockSite.  If we are in a block of solids, move on.
            int *lProcRankForSite =
                iGlobLatDat.Blocks[lBlockNumber].ProcessorRankForEachBlockSite;

            // If the array of proc rank for each site is NULL, we're on an all-solid block.
            if (lProcRankForSite == NULL)
            {
              continue;
            }

            // Create variables for the index of this site on the block and the number of fluid sites
            // that have been assigned to the current processor.
            int lSiteNumber = -1;

            // For each dimension of the site co-ordinates, iterate over all values of the site
            // co-ordinates on the current block.
            for (int lSiteCoordI = lBlockCoordI * iGlobLatDat.GetBlockSize(); lSiteCoordI
                < lBlockCoordI * iGlobLatDat.GetBlockSize()
                    + iGlobLatDat.GetBlockSize(); lSiteCoordI++)
            {
              for (int lSiteCoordJ = lBlockCoordJ * iGlobLatDat.GetBlockSize(); lSiteCoordJ
                  < lBlockCoordJ * iGlobLatDat.GetBlockSize()
                      + iGlobLatDat.GetBlockSize(); lSiteCoordJ++)
              {
                for (int lSiteCoordK = lBlockCoordK
                    * iGlobLatDat.GetBlockSize(); lSiteCoordK < lBlockCoordK
                    * iGlobLatDat.GetBlockSize() + iGlobLatDat.GetBlockSize(); lSiteCoordK++)
                {
                  // Keep track of the site number.
                  lSiteNumber++;

                  // Move on if the site is solid (ProcessorRankForEachBlockSite = BIG_NUMBER2) or has
                  // already been assigned to a rank (0 <= ProcessorRankForEachBlockSite < BIG_NUMBER2).
                  if (lProcRankForSite[lSiteNumber] != iMarker)
                  {
                    continue;
                  }
                  // We have found an unvisited fluid site to start growing the subdomain from.
                  // Assign it to the rank and update the fluid site counters.
                  lProcRankForSite[lSiteNumber] = proc_count;

                  ++lSitesOnCurrentProc;

                  // Record the location of this initial site.
                  lSiteLocationA->clear();
                  SiteLocation *lNew = new SiteLocation();
                  lNew->i = lSiteCoordI;
                  lNew->j = lSiteCoordJ;
                  lNew->k = lSiteCoordK;
                  lSiteLocationA->push_back(lNew);

                  // The subdomain can grow.
                  bool lIsRegionGrowing = true;

                  // While the region can grow (i.e. it is not bounded by solids or visited
                  // sites), and we need more sites on this particular rank.
                  while (lSitesOnCurrentProc < iSitesPerProc
                      && lIsRegionGrowing)
                  {
                    for (unsigned int ii = 0; ii < lSiteLocationB->size(); ii++)
                    {
                      delete lSiteLocationB->operator [](ii);
                    }
                    lSiteLocationB->clear();

                    // Sites added to the edge of the mClusters during the iteration.
                    lIsRegionGrowing = false;

                    // For sites on the edge of the domain (sites_a), deal with the neighbours.
                    for (unsigned int index_a = 0; index_a
                        < lSiteLocationA->size() && lSitesOnCurrentProc
                        < iSitesPerProc; index_a++)
                    {
                      lNew = lSiteLocationA->operator [](index_a);

                      for (unsigned int l = 1; l < D3Q15::NUMVECTORS
                          && lSitesOnCurrentProc < iSitesPerProc; l++)
                      {
                        // Record neighbour location.
                        int neigh_i = lNew->i + D3Q15::CX[l];
                        int neigh_j = lNew->j + D3Q15::CY[l];
                        int neigh_k = lNew->k + D3Q15::CZ[l];

                        // Move on if neighbour is outside the bounding box.
                        if (neigh_i == -1 || neigh_i
                            == iGlobLatDat.GetXSiteCount())
                          continue;
                        if (neigh_j == -1 || neigh_j
                            == iGlobLatDat.GetYSiteCount())
                          continue;
                        if (neigh_k == -1 || neigh_k
                            == iGlobLatDat.GetZSiteCount())
                          continue;

                        // Move on if the neighbour is in a block of solids (in which case
                        // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid or has already
                        // been assigned to a rank (in which case ProcessorRankForEachBlockSite != -1).  ProcessorRankForEachBlockSite
                        // was initialized in lbmReadConfig in io.cc.

                        // Pointer to the rank on which a particular fluid site
                        // resides.
                        int * proc_id_p =
                            iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i,
                                                                  neigh_j,
                                                                  neigh_k);

                        if (proc_id_p == NULL || *proc_id_p != iMarker)
                        {
                          continue;
                        }
                        // Set the rank for a neighbour and update the fluid site counters.
                        *proc_id_p = proc_count;
                        lSitesOnCurrentProc++;

                        // Neighbour was found, so the region can grow.
                        lIsRegionGrowing = true;

                        // Record the location of the neighbour.
                        SiteLocation * lNewB = new SiteLocation();
                        lNewB->i = neigh_i;
                        lNewB->j = neigh_j;
                        lNewB->k = neigh_k;
                        lSiteLocationB->push_back(lNewB);
                      }
                    }
                    // When the new layer of edge sites has been found, swap the buffers for
                    // the current and new layers of edge sites.
                    std::vector<SiteLocation*> *tempP = lSiteLocationA;
                    lSiteLocationA = lSiteLocationB;
                    lSiteLocationB = tempP;
                  }

                  // If we have enough sites, we have finished.
                  if (lSitesOnCurrentProc >= iSitesPerProc)
                  {
                    if (!iIsMachineLevel)
                    {
                      bNetTopology->FluidSitesOnEachProcessor[proc_count]
                          = lSitesOnCurrentProc;
                    }

                    ++proc_count;
                    if (!iIsMachineLevel)
                    {
                      bUnassignedSites -= lSitesOnCurrentProc;
                      iSitesPerProc = (int) ceil((double) bUnassignedSites
                          / (double) (bNetTopology->GetProcessorCount()
                              - proc_count));
                    }
                    lSitesOnCurrentProc = 0;
                  }
                  // If not, we have to start growing a different region for the same rank:
                  // region expansions could get trapped.

                } // Site co-ord k
              } // Site co-ord j
            } // Site co-ord i
          } // Block co-ord k
        } // Block co-ord j
      } // Block co-ord i

      for (unsigned int ii = 0; ii < lSiteLocationA->size(); ii++)
      {
        delete lSiteLocationA->operator [](ii);
      }
      for (unsigned int ii = 0; ii < lSiteLocationB->size(); ii++)
      {
        delete lSiteLocationB->operator [](ii);
      }

      delete lSiteLocationA;
      delete lSiteLocationB;
    }
  }
}
