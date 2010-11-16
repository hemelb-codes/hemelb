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
                                                           int & fluid_sites_per_unit,
                                                           int & unvisited_fluid_sites,
                                                           const int iCurrentProcId,
                                                           const int unitLevel,
                                                           lb::LocalLatticeData * iLocalLatDat,
                                                           lb::GlobalLatticeData &iGlobLatDat,
                                                           NetworkTopology * bNetTopology)
    {

      int sites_buffer_size = 10000;
      SiteLocation *site_location_a = new SiteLocation[sites_buffer_size];
      SiteLocation *site_location_b = new SiteLocation[sites_buffer_size];

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
            int lFluidSitesOnCurrentProcessor = 0;

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
                  if (lProcRankForSite[lSiteNumber] != iCurrentProcId)
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
                  ++lFluidSitesOnCurrentProcessor;

                  if (unitLevel == 0)
                  {
                    ++bNetTopology->FluidSitesOnEachProcessor[proc_count];
                  }

                  // Sites on the edge of the mClusters at the start of the current graph growing partitioning step.
                  int sites_a = 1;

                  // Record the location of this initial site.
                  SiteLocation *site_location_a_p = &site_location_a[0];
                  site_location_a_p->i = lSiteCoordI;
                  site_location_a_p->j = lSiteCoordJ;
                  site_location_a_p->k = lSiteCoordK;

                  // The subdomain can grow.
                  bool lAreFluidSitesIncrementing = true;

                  // While the region can grow (i.e. it is not bounded by solids or visited
                  // sites), and we need more sites on this particular rank.
                  while (lFluidSitesOnCurrentProcessor < fluid_sites_per_unit
                      && lAreFluidSitesIncrementing)
                  {
                    // Sites added to the edge of the mClusters during the iteration.
                    int sites_b = 0;
                    lAreFluidSitesIncrementing = false;

                    // For sites on the edge of the domain (sites_a), deal with the neighbours.
                    for (int index_a = 0; index_a < sites_a
                        && lFluidSitesOnCurrentProcessor < fluid_sites_per_unit; index_a++)
                    {
                      site_location_a_p = &site_location_a[index_a];

                      for (unsigned int l = 1; l < D3Q15::NUMVECTORS
                          && lFluidSitesOnCurrentProcessor
                              < fluid_sites_per_unit; l++)
                      {
                        // Record neighbour location.
                        int neigh_i = site_location_a_p->i + D3Q15::CX[l];
                        int neigh_j = site_location_a_p->j + D3Q15::CY[l];
                        int neigh_k = site_location_a_p->k + D3Q15::CZ[l];

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

                        if (proc_id_p == NULL || *proc_id_p != iCurrentProcId)
                        {
                          continue;
                        }
                        // Set the rank for a neighbour and update the fluid site counters.
                        *proc_id_p = proc_count;
                        lFluidSitesOnCurrentProcessor++;

                        if (unitLevel == 0)
                        {
                          bNetTopology->FluidSitesOnEachProcessor[proc_count]++;
                        }

                        if (bNetTopology->LocalRank == proc_count)
                        {
                          ++iLocalLatDat->LocalFluidSites;
                        }

                        // Neighbour was found, so the region can grow.
                        lAreFluidSitesIncrementing = true;

                        // If the new layer of neighbours is too large, allocate more
                        // memory.
                        if (sites_b == sites_buffer_size)
                        {
                          sites_buffer_size *= 2;
                          site_location_a
                              = (SiteLocation *) realloc(
                                                         site_location_a,
                                                         sizeof(SiteLocation)
                                                             * sites_buffer_size);
                          site_location_b
                              = (SiteLocation *) realloc(
                                                         site_location_b,
                                                         sizeof(SiteLocation)
                                                             * sites_buffer_size);
                        }

                        // Record the location of the neighbour.
                        SiteLocation * site_location_b_p =
                            &site_location_b[sites_b];
                        site_location_b_p->i = neigh_i;
                        site_location_b_p->j = neigh_j;
                        site_location_b_p->k = neigh_k;
                        ++sites_b;
                      }
                    }
                    // When the new layer of edge sites has been found, swap the buffers for
                    // the current and new layers of edge sites.
                    site_location_a_p = site_location_a;
                    site_location_a = site_location_b;
                    site_location_b = site_location_a_p;
                    sites_a = sites_b;
                  }

                  // If we have enough sites, we have finished.
                  if (lFluidSitesOnCurrentProcessor >= fluid_sites_per_unit)
                  {
                    ++proc_count;
                    if (unitLevel == 0)
                    {
                      unvisited_fluid_sites -= lFluidSitesOnCurrentProcessor;
                      fluid_sites_per_unit
                          = (int) ceil((double) unvisited_fluid_sites
                              / (double) (bNetTopology->ProcessorCount
                                  - proc_count));
                    }
                    lFluidSitesOnCurrentProcessor = 0;
                  }
                  // If not, we have to start growing a different region for the same rank:
                  // region expansions could get trapped.
                }
              }
            }
          }
        }
      }

      delete[] site_location_b;
      delete[] site_location_a;
    }
  }
}
