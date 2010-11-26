#include "topology/BaseTopologyManager.h"

namespace hemelb
{
  namespace topology
  {
    BaseTopologyManager::BaseTopologyManager()
    {
      // This exists to prevent instantiation.
    }

    void BaseTopologyManager::DecomposeDomain(int iTotalFluidSites,
                                              NetworkTopology & bNetTop,
                                              const lb::GlobalLatticeData & bGlobLatDat,
                                              lb::LocalLatticeData & bLocalLatDat)
    {
      // Allocations.  fluid sites will store actual number of fluid
      // sites per proc.  Site location will store up to 10000 of some
      // sort of coordinate.
      bLocalLatDat.LocalFluidSites = 0;

      // Initialise the count of fluid sites on each processor to 0.
      bNetTop.FluidSitesOnEachProcessor = new int[bNetTop.ProcessorCount];

      for (int n = 0; n < bNetTop.ProcessorCount; n++)
      {
        bNetTop.FluidSitesOnEachProcessor[n] = 0;
      }

      // Count of fluid sites not yet visited.
      int lUnvisitedFluidSiteCount = iTotalFluidSites;

      // If one machine or one machine per proc.
      if (bNetTop.MachineCount == 1 || bNetTop.MachineCount
          == bNetTop.ProcessorCount)
      {
        // Fluid sites per rank.
        int iSitesPerProc = (int) ceil((double) iTotalFluidSites
            / (double) bNetTop.ProcessorCount);

        //Rank we're looking at.
        int proc_count = 0;

        // If we're steering with more than one processor, save one processor for doing that.
#ifndef NO_STEER
        if (bNetTop.ProcessorCount != 1)
        {
          iSitesPerProc = (int) ceil((double) iTotalFluidSites
              / (double) (bNetTop.ProcessorCount - 1));
          proc_count = 1;
        }
#endif

        // In the simple case, simply divide fluid sites up between processors.
        AssignFluidSitesToProcessors(proc_count, iSitesPerProc,
                                     lUnvisitedFluidSiteCount, -1, false,
                                     &bLocalLatDat, bGlobLatDat, &bNetTop);
      }
      else
      {
        // Rank we are looking at.
        int proc_count = bNetTop.ProcessorCount;
        double weight = (double) (bNetTop.ProcCountOnEachMachine[0]
            * bNetTop.ProcessorCount) / (double) (bNetTop.ProcessorCount - 1);
        // Fluid sites per rank.
        int iSitesPerProc = (int) ceil((double) iTotalFluidSites * weight
            / bNetTop.MachineCount);

        // First, divide the sites up between machines.
        AssignFluidSitesToProcessors(proc_count, iSitesPerProc,
                                     lUnvisitedFluidSiteCount, -1, true,
                                     &bLocalLatDat, bGlobLatDat, &bNetTop);

        iSitesPerProc = (int) ceil((double) lUnvisitedFluidSiteCount
            / (double) (bNetTop.ProcessorCount - 1));
        proc_count = 1;

        // For each machine, divide up the sites it has between its cores.
        for (int lMachineNumber = 0; lMachineNumber < bNetTop.MachineCount; lMachineNumber++)
        {
          AssignFluidSitesToProcessors(proc_count, iSitesPerProc,
                                       lUnvisitedFluidSiteCount,
                                       bNetTop.ProcessorCount + lMachineNumber,
                                       false, &bLocalLatDat, bGlobLatDat,
                                       &bNetTop);
        }
      }

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
      int partial_visited_fluid_sites = 0;

      int n = -1;
      hemelb::lb::BlockData* proc_block_p;

      int site_i, site_j, site_k; // Global coordinates of a site.
      int neigh_i, neigh_j, neigh_k; // Global coordinates of a neighbour site.
      int i, j, k; // Global block index.
      int l; // Index for neighbours of a site.
      int m; // Site index on a paricular block.
      int mm; // Index of processors surrounding this one.
      int are_fluid_sites_incrementing;
      int block_size = iGlobLatDat.BlockSize;
      int sites_a, sites_b, index_a;
      int sites_buffer_size = 10000;
      SiteLocation *site_location_a, *site_location_b;
      SiteLocation *site_location_a_p, *site_location_b_p;
      site_location_a = (SiteLocation *) malloc(sizeof(SiteLocation)
          * sites_buffer_size);
      site_location_b = (SiteLocation *) malloc(sizeof(SiteLocation)
          * sites_buffer_size);

      int sites_x = iGlobLatDat.SitesX;
      int sites_y = iGlobLatDat.SitesY;
      int sites_z = iGlobLatDat.SitesZ;

      int* proc_id_p;

      for (i = 0; i < iGlobLatDat.BlocksX; i++)
        for (j = 0; j < iGlobLatDat.BlocksY; j++)
          for (k = 0; k < iGlobLatDat.BlocksZ; k++)
          {
            proc_block_p = &iGlobLatDat.Blocks[++n];

            if (proc_block_p->ProcessorRankForEachBlockSite == NULL)
            {
              continue;
            }
            m = -1;

            for (site_i = i * block_size; site_i < i * block_size + block_size; site_i++)
              for (site_j = j * block_size; site_j < j * block_size
                  + block_size; site_j++)
                for (site_k = k * block_size; site_k < k * block_size
                    + block_size; site_k++)
                {
                  if (proc_block_p->ProcessorRankForEachBlockSite[++m]
                      != iMarker)
                  {
                    continue;
                  }
                  proc_block_p->ProcessorRankForEachBlockSite[m] = proc_count;

                  if (proc_count == bNetTopology->LocalRank)
                  {
                    ++iLocalLatDat->LocalFluidSites;
                  }
                  ++partial_visited_fluid_sites;

                  if (!iIsMachineLevel)
                    ++bNetTopology->FluidSitesOnEachProcessor[proc_count];

                  sites_a = 1;
                  site_location_a_p = &site_location_a[0];
                  site_location_a_p->i = site_i;
                  site_location_a_p->j = site_j;
                  site_location_a_p->k = site_k;

                  are_fluid_sites_incrementing = 1;

                  while (partial_visited_fluid_sites < iSitesPerProc
                      && are_fluid_sites_incrementing)
                  {
                    sites_b = 0;
                    are_fluid_sites_incrementing = 0;

                    for (index_a = 0; index_a < sites_a
                        && partial_visited_fluid_sites < iSitesPerProc; index_a++)
                    {
                      site_location_a_p = &site_location_a[index_a];

                      for (l = 1; l < 15 && partial_visited_fluid_sites
                          < iSitesPerProc; l++)
                      {
                        neigh_i = site_location_a_p->i + D3Q15::CX[l];
                        neigh_j = site_location_a_p->j + D3Q15::CY[l];
                        neigh_k = site_location_a_p->k + D3Q15::CZ[l];

                        if (neigh_i == -1 || neigh_i == sites_x)
                          continue;
                        if (neigh_j == -1 || neigh_j == sites_y)
                          continue;
                        if (neigh_k == -1 || neigh_k == sites_z)
                          continue;

                        proc_id_p
                            = iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i,
                                                                    neigh_j,
                                                                    neigh_k);

                        if (proc_id_p == NULL || *proc_id_p != iMarker)
                        {
                          continue;
                        }
                        *proc_id_p = proc_count;

                        ++partial_visited_fluid_sites;

                        if (!iIsMachineLevel)
                          ++bNetTopology->FluidSitesOnEachProcessor[proc_count];

                        are_fluid_sites_incrementing = 1;

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
                        site_location_b_p = &site_location_b[sites_b];
                        site_location_b_p->i = neigh_i;
                        site_location_b_p->j = neigh_j;
                        site_location_b_p->k = neigh_k;
                        ++sites_b;

                        if (proc_count == bNetTopology->LocalRank)
                        {
                          ++iLocalLatDat->LocalFluidSites;
                        }
                      }
                    }
                    site_location_a_p = site_location_a;
                    site_location_a = site_location_b;
                    site_location_b = site_location_a_p;
                    sites_a = sites_b;
                  }
                  if (partial_visited_fluid_sites >= iSitesPerProc)
                  {
                    ++proc_count;

                    if (!iIsMachineLevel)
                    {
                      bUnassignedSites -= partial_visited_fluid_sites;
                      iSitesPerProc
                          = (int) ceil((double) bUnassignedSites
                              / (double) (bNetTopology->ProcessorCount
                                  - proc_count));
                    }
                    partial_visited_fluid_sites = 0;
                  }
                }
          }
      free(site_location_b);
      free(site_location_a);
    }
  }
}
