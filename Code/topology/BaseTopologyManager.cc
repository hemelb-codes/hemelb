#include "topology/BaseTopologyManager.h"

namespace hemelb
{
  namespace topology
  {
    BaseTopologyManager::BaseTopologyManager()
    {
      // This exists to prevent instantiation.
    }

    void BaseTopologyManager::AssignFluidSitesToProcessors(int & proc_count,
                                                           int & iSitesPerProc,
                                                           int & bUnassignedSites,
                                                           const int iMarker,
                                                           const bool iIsMachineLevel,
                                                           lb::LocalLatticeData * iLocalLatDat,
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

      // Iterate over all blocks.
      for (int lBlockCoordI = 0; lBlockCoordI < iGlobLatDat.BlocksX; lBlockCoordI++)
      {
        for (int lBlockCoordJ = 0; lBlockCoordJ < iGlobLatDat.BlocksY; lBlockCoordJ++)
        {
          for (int lBlockCoordK = 0; lBlockCoordK < iGlobLatDat.BlocksZ; lBlockCoordK++)
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
            int lSitesOnCurrentProc = 0;

            // For each dimension of the site co-ordinates, iterate over all values of the site
            // co-ordinates on the current block.
            for (int lSiteCoordI = lBlockCoordI * iGlobLatDat.BlockSize; lSiteCoordI
                < lBlockCoordI * iGlobLatDat.BlockSize + iGlobLatDat.BlockSize; lSiteCoordI++)
            {
              for (int lSiteCoordJ = lBlockCoordJ * iGlobLatDat.BlockSize; lSiteCoordJ
                  < lBlockCoordJ * iGlobLatDat.BlockSize
                      + iGlobLatDat.BlockSize; lSiteCoordJ++)
              {
                for (int lSiteCoordK = lBlockCoordK * iGlobLatDat.BlockSize; lSiteCoordK
                    < lBlockCoordK * iGlobLatDat.BlockSize
                        + iGlobLatDat.BlockSize; lSiteCoordK++)
                {
                  // Keep track of the site number.
                  lSiteNumber++;

                  //TODO comments from here.
                  // Move on if the site is solid (ProcessorRankForEachBlockSite = 1 << 30) or has
                  // already been assigned to a rank (0 <= ProcessorRankForEachBlockSite < 1 << 30).
                  if (lProcRankForSite[lSiteNumber] != iMarker)
                  {
                    continue;
                  }
                  // We have found an unvisited fluid site to start growing the subdomain from.
                  // Assign it to the rank and update the fluid site counters.
                  lProcRankForSite[lSiteNumber] = proc_count;

                  if (bNetTopology->LocalRank == proc_count)
                  {
                    ++iLocalLatDat->LocalFluidSites;
                  }
                  ++lSitesOnCurrentProc;

                  if (!iIsMachineLevel)
                  {
                    ++bNetTopology->FluidSitesOnEachProcessor[proc_count];
                  }

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
                    for(unsigned int ii = 0; ii < lSiteLocationB->size(); ii++)
                    {
                      delete lSiteLocationB->operator [](ii);
                    }
                    lSiteLocationB->clear();

                    // Sites added to the edge of the mClusters during the iteration.
                    lIsRegionGrowing = false;

                    // For sites on the edge of the domain (sites_a), deal with the neighbours.
                    for (unsigned int index_a = 0; index_a < lSiteLocationA->size()
                        && lSitesOnCurrentProc < iSitesPerProc; index_a++)
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
                        if (neigh_i == -1 || neigh_i == iGlobLatDat.SitesX)
                          continue;
                        if (neigh_j == -1 || neigh_j == iGlobLatDat.SitesY)
                          continue;
                        if (neigh_k == -1 || neigh_k == iGlobLatDat.SitesZ)
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

                        if (!iIsMachineLevel)
                        {
                          bNetTopology->FluidSitesOnEachProcessor[proc_count]++;
                        }

                        if (bNetTopology->LocalRank == proc_count)
                        {
                          ++iLocalLatDat->LocalFluidSites;
                        }

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
                    ++proc_count;
                    if (!iIsMachineLevel)
                    {
                      bUnassignedSites -= lSitesOnCurrentProc;
                      iSitesPerProc
                          = (int) ceil((double) bUnassignedSites
                              / (double) (bNetTopology->ProcessorCount
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

      for(unsigned int ii = 0; ii < lSiteLocationA->size(); ii++)
      {
        delete lSiteLocationA->operator [](ii);
      }
      for(unsigned int ii = 0; ii < lSiteLocationB->size(); ii++)
      {
        delete lSiteLocationB->operator [](ii);
      }

      delete lSiteLocationA;
      delete lSiteLocationB;
    }
  }
}
