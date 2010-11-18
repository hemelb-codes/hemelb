/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include "lb.h"
#include "net.h"
#include "utilityFunctions.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

// Returns the type of collision/streaming update for the fluid site
// with data "site_data".
unsigned int Net::GetCollisionType(unsigned int site_data)
{
  unsigned int boundary_type;

  if (site_data == hemelb::lb::FLUID_TYPE)
  {
    return FLUID;
  }
  boundary_type = site_data & SITE_TYPE_MASK;

  if (boundary_type == hemelb::lb::FLUID_TYPE)
  {
    return EDGE;
  }
  if (! (site_data & PRESSURE_EDGE_MASK))
  {
    if (boundary_type == hemelb::lb::INLET_TYPE)
    {
      return INLET;
    }
    else
    {
      return OUTLET;
    }
  }
  else
  {
    if (boundary_type == hemelb::lb::INLET_TYPE)
    {
      return INLET | EDGE;
    }
    else
    {
      return OUTLET | EDGE;
    }
  }
}

void Net::Abort()
{
#ifndef NOMPI
  err = MPI_Abort(MPI_COMM_WORLD, 1);
#else
  exit(1);
#endif
}

/*!
 This is called from the main function.  First function to deal with processors.
 The domain partitioning technique and the management of the
 buffers useful for the inter-processor communications are
 implemented in this function.  The domain decomposition is based
 on a graph growing partitioning technique.
 */
void Net::Initialise(int iTotalFluidSites,
                     hemelb::topology::TopologyManager &iTopologyManager,
                     hemelb::topology::NetworkTopology &iNetTop,
                     hemelb::lb::GlobalLatticeData &iGlobLatDat,
                     hemelb::lb::LocalLatticeData &bLocalLatDat)
{
  block_count = iGlobLatDat.BlockCount;

  // Allocations.  fluid sites will store actual number of fluid
  // sites per proc.  Site location will store up to 10000 of some
  // sort of coordinate.
  iNetTop.FluidSitesOnEachProcessor = new int[iNetTop.ProcessorCount];
  bLocalLatDat.LocalFluidSites = 0;

  // a fast graph growing partitioning technique which spans the data
  // set only once is implemented here; the data set is explored by
  // means of the arrays "site_location_a[]" and "site_location_b[]"

  for (int n = 0; n < iNetTop.ProcessorCount; n++)
  {
    iNetTop.FluidSitesOnEachProcessor[n] = 0;
  }

  // Count of fluid sites not yet visited.
  int lUnvisitedFluidSiteCount = iTotalFluidSites;

  // Count of time elapsed.
  double seconds = hemelb::util::myClock();

  // If one machine or one machine per proc.
  if (iNetTop.MachineCount == 1 || iNetTop.MachineCount
      == iNetTop.ProcessorCount)
  {
    // Fluid sites per rank.
    int fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
        / (double) iNetTop.ProcessorCount);

    //Rank we're looking at.
    int proc_count = 0;

    // If we're steering with more than one processor, save one processor for doing that.
#ifndef NO_STEER
    if (iNetTop.ProcessorCount != 1)
    {
      fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
          / (double) (iNetTop.ProcessorCount - 1));
      proc_count = 1;
    }
#endif

    // In the simple case, simply divide fluid sites up between processors.
    iTopologyManager.AssignFluidSitesToProcessors(proc_count,
                                                  fluid_sites_per_unit,
                                                  lUnvisitedFluidSiteCount, -1,
                                                  false, &bLocalLatDat,
                                                  iGlobLatDat, &iNetTop);
  }
  else
  {
    // Rank we are looking at.
    int proc_count = iNetTop.ProcessorCount;
    double weight = (double) (iNetTop.ProcCountOnEachMachine[0]
        * iNetTop.ProcessorCount) / (double) (iNetTop.ProcessorCount - 1);
    // Fluid sites per rank.
    int fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites * weight
        / iNetTop.MachineCount);

    // First, divide the sites up between machines.
    iTopologyManager.AssignFluidSitesToProcessors(proc_count,
                                                  fluid_sites_per_unit,
                                                  lUnvisitedFluidSiteCount, -1,
                                                  true, &bLocalLatDat,
                                                  iGlobLatDat, &iNetTop);

    fluid_sites_per_unit = (int) ceil((double) lUnvisitedFluidSiteCount
        / (double) (iNetTop.ProcessorCount - 1));
    proc_count = 1;

    // For each machine, divide up the sites it has between its cores.
    for (int lMachineNumber = 0; lMachineNumber < iNetTop.MachineCount; lMachineNumber++)
    {
      iTopologyManager.AssignFluidSitesToProcessors(proc_count,
                                                    fluid_sites_per_unit,
                                                    lUnvisitedFluidSiteCount,
                                                    iNetTop.ProcessorCount
                                                        + lMachineNumber, false,
                                                    &bLocalLatDat, iGlobLatDat,
                                                    &iNetTop);
    }
  }

  // Next stage of the timing.
  dd_time = hemelb::util::myClock() - seconds;
  seconds = hemelb::util::myClock();

  // Create a map between the two-level data representation and the 1D
  // compact one is created here.

  // This rank's site data.
  unsigned int *lThisRankSiteData =
      new unsigned int[bLocalLatDat.LocalFluidSites];

  // Array of booleans to store whether any sites on a block are fluid
  // sites residing on this rank.
  bool *lBlockIsOnThisRank = new bool[iGlobLatDat.BlockCount];
  // Initialise to false.
  for (int n = 0; n < iGlobLatDat.BlockCount; n++)
  {
    lBlockIsOnThisRank[n] = false;
  }

  int lSiteIndexOnProc = 0;

  for (int lBlockNumber = 0; lBlockNumber < iGlobLatDat.BlockCount; lBlockNumber++)
  {
    hemelb::lb::BlockData * lCurrentDataBlock =
        &iGlobLatDat.Blocks[lBlockNumber];

    // If we are in a block of solids, move to the next block.
    if (lCurrentDataBlock->site_data == NULL)
    {
      continue;
    }

    // If we have some fluid sites, point to mProcessorsForEachBlock and map_block.
    hemelb::lb::BlockData * proc_block_p = &iGlobLatDat.Blocks[lBlockNumber];

    // lCurrentDataBlock.site_data is set to the fluid site identifier on this rank or (1U << 31U) if a site is solid
    // or not on this rank.  site_data is indexed by fluid site identifier and set to the site_data.
    for (int lSiteIndexWithinBlock = 0; lSiteIndexWithinBlock
        < iGlobLatDat.SitesPerBlockVolumeUnit; lSiteIndexWithinBlock++)
    {
      if (iNetTop.LocalRank
          == proc_block_p->ProcessorRankForEachBlockSite[lSiteIndexWithinBlock])
      {
        // If the current site is non-solid, copy the site data into the array for
        // this rank (in the whole-processor location), then set the site data
        // for this site within the current block to be the site index over the whole
        // processor.
        if ( (lCurrentDataBlock->site_data[lSiteIndexWithinBlock]
            & SITE_TYPE_MASK) != hemelb::lb::SOLID_TYPE)
        {
          lThisRankSiteData[lSiteIndexOnProc]
              = lCurrentDataBlock->site_data[lSiteIndexWithinBlock];
          lCurrentDataBlock->site_data[lSiteIndexWithinBlock]
              = lSiteIndexOnProc;
          ++lSiteIndexOnProc;
        }
        else
        {
          // If this is a solid, set the site data on the current block to
          // some massive value.
          lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = (1U << 31U);
        }
        // Set the array to notify that the current block has sites on this
        // rank.
        lBlockIsOnThisRank[lBlockNumber] = true;
      }
      // If this site is not on the current processor, set its whole processor
      // index within the per-block store to a nonsense value.
      else
      {
        lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = (1U << 31U);
      }
    }
  }

  // If we are in a block of solids, we set map_block[n].site_data to NULL.
  for (int n = 0; n < iGlobLatDat.BlockCount; n++)
  {
    if (lBlockIsOnThisRank[n])
    {
      continue;
    }

    delete[] iGlobLatDat.Blocks[n].site_data;
    iGlobLatDat.Blocks[n].site_data = NULL;

    if (iGlobLatDat.Blocks[n].wall_data != NULL)
    {
      delete[] iGlobLatDat.Blocks[n].wall_data;
      iGlobLatDat.Blocks[n].wall_data = NULL;
    }
  }
  delete[] lBlockIsOnThisRank;

  // The numbers of inter- and intra-machine neighbouring processors,
  // interface-dependent and independent fluid sites and shared
  // distribution functions of the reference processor are calculated
  // here.  neigh_proc is a static array that is declared in config.h.

  // Initialise various things to 0.
  my_inter_sites = 0;
  my_inner_sites = 0;

  for (int m = 0; m < COLLISION_TYPES; m++)
  {
    my_inter_collisions[m] = 0;
    my_inner_collisions[m] = 0;
  }

  iNetTop.TotalSharedFs = 0; // shared SharedFCount within Net struct.

  lSiteIndexOnProc = 0;

  int n = -1;

  // Iterate over all blocks in site units
  for (int i = 0; i < iGlobLatDat.SitesX; i += iGlobLatDat.BlockSize)
  {
    for (int j = 0; j < iGlobLatDat.SitesY; j += iGlobLatDat.BlockSize)
    {
      for (int k = 0; k < iGlobLatDat.SitesZ; k += iGlobLatDat.BlockSize)
      {
        hemelb::lb::BlockData * map_block_p = &iGlobLatDat.Blocks[++n];

        if (map_block_p->site_data == NULL)
        {
          continue;
        }

        int m = -1;

        // Iterate over all sites within the current block.
        for (int site_i = i; site_i < i + iGlobLatDat.BlockSize; site_i++)
        {
          for (int site_j = j; site_j < j + iGlobLatDat.BlockSize; site_j++)
          {
            for (int site_k = k; site_k < k + iGlobLatDat.BlockSize; site_k++)
            {
              m++;
              // If the site is not on this processor, continue.
              if (iNetTop.LocalRank
                  != map_block_p->ProcessorRankForEachBlockSite[m])
              {
                continue;
              }

              bool lIsInnerSite = true;

              // Iterate over all direction vectors.
              for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
              {
                // Find the neighbour site co-ords in this direction.
                int neigh_i = site_i + D3Q15::CX[l];
                int neigh_j = site_j + D3Q15::CY[l];
                int neigh_k = site_k + D3Q15::CZ[l];

                // Find the processor Id for that neighbour.
                int *proc_id_p = iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i,
                                                                       neigh_j,
                                                                       neigh_k);

                // Move on if the neighbour is in a block of solids (in which case
                // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                // 1 << 30) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                // in lbmReadConfig in io.cc.
                if (proc_id_p == NULL || iNetTop.LocalRank == (*proc_id_p)
                    || *proc_id_p == (1 << 30))
                {
                  continue;
                }

                lIsInnerSite = false;

                // The first time, we set mm = 0, flag = 1, but net_neigh_procs = 0, so
                // the loop is not executed.
                int mm, flag;

                // Iterate over neighbouring processors until we find the one with the
                // neighbouring site on it.
                int lNeighbouringProcs = iNetTop.NeighbouringProcs.size();
                for (mm = 0, flag = 1; mm < lNeighbouringProcs && flag; mm++)
                {
                  // Check whether the rank for a particular neighbour has already been
                  // used for this processor.  If it has, set flag to zero.
                  hemelb::topology::NeighbouringProcessor * neigh_proc_p =
                      iNetTop.NeighbouringProcs[mm];

                  // If ProcessorRankForEachBlockSite is equal to a neigh_proc that has alredy been listed.
                  if (*proc_id_p == neigh_proc_p->Rank)
                  {
                    flag = 0;
                    ++neigh_proc_p->SharedFCount;
                    ++iNetTop.TotalSharedFs;
                  }
                }
                // If flag is 1, we didn't find a neighbour-proc with the neighbour-site on it
                // so we need a new neighbouring processor.
                if (flag)
                {
                  // Store rank of neighbour in >neigh_proc[neigh_procs]
                  hemelb::topology::NeighbouringProcessor * lNewNeighbour =
                      new hemelb::topology::NeighbouringProcessor();
                  lNewNeighbour->SharedFCount = 1;
                  lNewNeighbour->Rank = *proc_id_p;
                  iNetTop.NeighbouringProcs.push_back(lNewNeighbour);
                  ++iNetTop.TotalSharedFs;
                }
              }

              // Set the collision type data. map_block site data is renumbered according to
              // fluid site numbers within a particular collision type.

              int l = -1;

              if (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc])
                  == FLUID)
              {
                l = 0;
              }
              else if (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc])
                  == EDGE)
              {
                l = 1;
              }
              else if (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc])
                  == INLET)
              {
                l = 2;
              }
              else if (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc])
                  == OUTLET)
              {
                l = 3;
              }
              else if (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc])
                  == (INLET | EDGE))
              {
                l = 4;
              }
              else if (GetCollisionType(lThisRankSiteData[lSiteIndexOnProc])
                  == (OUTLET | EDGE))
              {
                l = 5;
              }

              ++lSiteIndexOnProc;

              if (lIsInnerSite)
              {
                ++my_inner_sites;

                if (l == 0)
                {
                  map_block_p->site_data[m] = my_inner_collisions[l];
                }
                else
                {
                  map_block_p->site_data[m] = 50000000 * (10 + (l - 1))
                      + my_inner_collisions[l];
                }
                ++my_inner_collisions[l];
              }
              else
              {
                ++my_inter_sites;

                if (l == 0)
                {
                  map_block_p->site_data[m] = 1000000000
                      + my_inter_collisions[l];
                }
                else
                {
                  map_block_p->site_data[m] = 50000000 * (20 + l)
                      + my_inter_collisions[l];
                }
                ++my_inter_collisions[l];
              }
            }
          }
        }
      }
    }
  }

  int collision_offset[2][COLLISION_TYPES];
  // Calculate the number of each type of collision.
  collision_offset[0][0] = 0;

  for (unsigned int l = 1; l < COLLISION_TYPES; l++)
  {
    collision_offset[0][l] = collision_offset[0][l - 1] + my_inner_collisions[l
        - 1];
  }
  collision_offset[1][0] = my_inner_sites;
  for (unsigned int l = 1; l < COLLISION_TYPES; l++)
  {
    collision_offset[1][l] = collision_offset[1][l - 1] + my_inter_collisions[l
        - 1];
  }

  // Iterate over blocks
  for (int n = 0; n < iGlobLatDat.BlockCount; n++)
  {
    hemelb::lb::BlockData *map_block_p = &iGlobLatDat.Blocks[n];

    // If we are in a block of solids, continue.
    if (map_block_p->site_data == NULL)
    {
      continue;
    }

    // Iterate over sites within the block.
    for (int m = 0; m < iGlobLatDat.SitesPerBlockVolumeUnit; m++)
    {
      unsigned int *site_data_p = &map_block_p->site_data[m];

      // If the site is solid, continue.
      if (*site_data_p & (1U << 31U))
      {
        continue;
      }

      // 0th collision type for inner sites, so don't do anything.
      if (*site_data_p < 500000000)
      {
        continue;
      }

      // Renumber the sites in map_block so that the numbers are compacted together.  We have
      // collision offset to tell us when one collision type ends and another starts.
      for (unsigned int l = 1; l < COLLISION_TYPES; l++)
      {
        if (*site_data_p >= 50000000 * (10 + (l - 1)) && *site_data_p
            < 50000000 * (10 + l))
        {
          *site_data_p += collision_offset[0][l] - 50000000 * (10 + (l - 1));
          break;
        }
      }
      for (unsigned int l = 0; l < COLLISION_TYPES; l++)
      {
        if (*site_data_p >= 50000000 * (20 + l) && *site_data_p < 50000000
            * (20 + (l + 1)))
        {
          *site_data_p += collision_offset[1][l] - 50000000 * (20 + l);
          break;
        }
      }
    }
  }

  // the precise interface-dependent data (interface-dependent fluid
  // site locations and identifiers of the distribution functions
  // streamed between different partitions) are collected and the
  // buffers needed for the communications are set from here

  short int *f_data = new short int[4 * iNetTop.TotalSharedFs];

  // Allocate the index in which to put the distribution functions received from the other
  // process.
  f_recv_iv = new int[iNetTop.TotalSharedFs];

  // Reset to zero again.
  iNetTop.TotalSharedFs = 0;

  // Set the remaining neighbouring processor data.
  for (unsigned int n = 0; n < iNetTop.NeighbouringProcs.size(); n++)
  {
    // f_data compacted according to number of shared f_s on each process.
    // f_data will be set later.
    iNetTop.NeighbouringProcs[n]->SharedFLocation
        = &f_data[iNetTop.TotalSharedFs << 2];

    // Pointing to a few things, but not setting any variables.
    // FirstSharedF points to start of shared_fs.
    iNetTop.NeighbouringProcs[n]->FirstSharedF = bLocalLatDat.LocalFluidSites
        * D3Q15::NUMVECTORS + 1 + iNetTop.TotalSharedFs;

    iNetTop.NeighbouringProcs[n]->SharedFReceivingIndex
        = &f_recv_iv[iNetTop.TotalSharedFs];

    iNetTop.TotalSharedFs += iNetTop.NeighbouringProcs[n]->SharedFCount;
    iNetTop.NeighbouringProcs[n]->SharedFCount = 0;// This is set back to 0.
  }

  iNetTop.NeighbourIndexFromProcRank = new short int[iNetTop.ProcessorCount];

  for (int m = 0; m < iNetTop.ProcessorCount; m++)
  {
    iNetTop.NeighbourIndexFromProcRank[m] = -1;
  }
  // Get neigh_proc_index from ProcessorRankForEachBlockSite.
  for (unsigned int m = 0; m < iNetTop.NeighbouringProcs.size(); m++)
  {
    iNetTop.NeighbourIndexFromProcRank[iNetTop.NeighbouringProcs[m]->Rank] = m;
  }

  // Allocate f_old and f_new according to the number of sites on the process.  The extra site
  // is there for when we would stream into a solid site during the simulation, which avoids
  // an if condition at every timestep at every boundary site.  We also allocate space for the
  // shared distribution functions.  We need twice as much space when we check the convergence
  // and the extra distribution functions are
  bLocalLatDat.FOld = new double[bLocalLatDat.LocalFluidSites
      * D3Q15::NUMVECTORS + 1 + iNetTop.TotalSharedFs];
  bLocalLatDat.FNew = new double[bLocalLatDat.LocalFluidSites
      * D3Q15::NUMVECTORS + 1 + iNetTop.TotalSharedFs];

  if (bLocalLatDat.LocalFluidSites > 0)
  {
    // f_id is allocated so we know which sites to get information from.
    bLocalLatDat.mFNeighbours = new int[bLocalLatDat.LocalFluidSites
        * D3Q15::NUMVECTORS];

    bLocalLatDat.mSiteData = new unsigned int[bLocalLatDat.LocalFluidSites];
    bLocalLatDat.mWallNormalAtSite = new double[bLocalLatDat.LocalFluidSites
        * 3];
    bLocalLatDat.mDistanceToWall = new double[bLocalLatDat.LocalFluidSites
        * (D3Q15::NUMVECTORS - 1)];
  }

  lSiteIndexOnProc = 0;

  n = -1;

  // Iterate over blocks in global co-ords.
  for (int i = 0; i < iGlobLatDat.SitesX; i += iGlobLatDat.BlockSize)
  {
    for (int j = 0; j < iGlobLatDat.SitesY; j += iGlobLatDat.BlockSize)
    {
      for (int k = 0; k < iGlobLatDat.SitesZ; k += iGlobLatDat.BlockSize)
      {
        hemelb::lb::BlockData *map_block_p = &iGlobLatDat.Blocks[++n];

        if (map_block_p->site_data == NULL)
        {
          continue;
        }

        hemelb::lb::BlockData *proc_block_p = &iGlobLatDat.Blocks[n];

        int m = -1;

        // Iterate over sites within the block.
        for (int site_i = i; site_i < i + iGlobLatDat.BlockSize; site_i++)
        {
          for (int site_j = j; site_j < j + iGlobLatDat.BlockSize; site_j++)
          {
            for (int site_k = k; site_k < k + iGlobLatDat.BlockSize; site_k++)
            {
              // If a site is not on this process, continue.
              m++;
              if (iNetTop.LocalRank
                  != proc_block_p->ProcessorRankForEachBlockSite[m])
              {
                continue;
              }
              // Get site data, which is the number of the fluid site on this proc..
              unsigned int site_map = map_block_p->site_data[m];

              // set f_id.
              bLocalLatDat.mFNeighbours[site_map * D3Q15::NUMVECTORS + 0]
                  = site_map * D3Q15::NUMVECTORS + 0;

              for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
              {
                // Work out positions of neighbours.
                int neigh_i = site_i + D3Q15::CX[l];
                int neigh_j = site_j + D3Q15::CY[l];
                int neigh_k = site_k + D3Q15::CZ[l];

                // initialize f_id to the rubbish site.
                bLocalLatDat.mFNeighbours[site_map * D3Q15::NUMVECTORS + l]
                    = bLocalLatDat.LocalFluidSites * D3Q15::NUMVECTORS;

                // You know which process the neighbour is on.
                int *proc_id_p = iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i,
                                                                       neigh_j,
                                                                       neigh_k);

                if (proc_id_p == NULL || *proc_id_p == 1 << 30)
                {
                  continue;
                }
                // Pointer to the neihgbour.
                unsigned int *site_data_p = iGlobLatDat.GetSiteData(neigh_i,
                                                                    neigh_j,
                                                                    neigh_k);

                // If on the same proc, set f_id of the
                // current site and direction to the
                // site and direction that it sends to.
                // If we check convergence, the data for
                // each site is split into that for the
                // current and previous cycles.
                if (iNetTop.LocalRank == *proc_id_p)
                {
                  bLocalLatDat.mFNeighbours[site_map * D3Q15::NUMVECTORS + l]
                      = *site_data_p * D3Q15::NUMVECTORS + l;

                  continue;
                }
                short int neigh_proc_index =
                    iNetTop.NeighbourIndexFromProcRank[*proc_id_p];

                // You have neigh proc again.
                hemelb::topology::NeighbouringProcessor * neigh_proc_p =
                    iNetTop.NeighbouringProcs[neigh_proc_index];

                // This stores some coordinates.  We
                // still need to know the site number.
                // neigh_proc[ n ].f_data is now
                // set as well, since this points to
                // f_data.  Every process has data for
                // its neighbours which say which sites
                // on this process are shared with the
                // neighbour.
                short int *f_data_p =
                    &neigh_proc_p->SharedFLocation[neigh_proc_p->SharedFCount
                        << 2];
                f_data_p[0] = site_i;
                f_data_p[1] = site_j;
                f_data_p[2] = site_k;
                f_data_p[3] = l;
                ++neigh_proc_p->SharedFCount; // We recount this again.
              }
              // This is used in Calculate BC in IO.
              bLocalLatDat.mSiteData[site_map]
                  = lThisRankSiteData[lSiteIndexOnProc];

              if (GetCollisionType(bLocalLatDat.mSiteData[site_map]) & EDGE)
              {
                for (unsigned int l = 0; l < 3; l++)
                  bLocalLatDat.mWallNormalAtSite[site_map * 3 + l]
                      = iGlobLatDat.Blocks[n].wall_data[m].wall_nor[l];

                for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
                  bLocalLatDat.mDistanceToWall[site_map * (D3Q15::NUMVECTORS
                      - 1) + l]
                      = iGlobLatDat.Blocks[n].wall_data[m].cut_dist[l];
              }
              else
              {
                bLocalLatDat.mWallNormalAtSite[site_map * 3] = BIG_NUMBER;
              }
              ++lSiteIndexOnProc;
            }
          }
        }
      }
    }
  }
  delete[] lThisRankSiteData;

  // point-to-point communications are performed to match data to be
  // sent to/receive from different partitions; in this way, the
  // communication of the locations of the interface-dependent fluid
  // sites and the identifiers of the distribution functions which
  // propagate to different partitions is avoided (only their values
  // will be communicated). It's here!

  // Allocate the request variable.
#ifndef NOMPI
  req = new MPI_Request*[COMMS_LEVELS];

  for (int m = 0; m < COMMS_LEVELS; m++)
  {
    req[m] = new MPI_Request[2 * iNetTop.ProcessorCount];
  }
#endif

  for (unsigned int m = 0; m < iNetTop.NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor * neigh_proc_p =
        iNetTop.NeighbouringProcs[m];

    // One way send receive.  The lower numbered iNetTop.ProcessorCount send and the higher numbered ones receive.
    // It seems that, for each pair of processors, the lower numbered one ends up with its own
    // edge sites and directions stored and the higher numbered one ends up with those on the
    // other processor.
#ifndef NOMPI
    if (neigh_proc_p->Rank > iNetTop.LocalRank)
    {
      err = MPI_Isend(&neigh_proc_p->SharedFLocation[0],
                      neigh_proc_p->SharedFCount * 4, MPI_SHORT,
                      neigh_proc_p->Rank, 10, MPI_COMM_WORLD, &req[0][m]);
    }
    else
    {
      err = MPI_Irecv(&neigh_proc_p->SharedFLocation[0],
                      neigh_proc_p->SharedFCount * 4, MPI_SHORT,
                      neigh_proc_p->Rank, 10, MPI_COMM_WORLD,
                      &req[0][iNetTop.NeighbouringProcs.size() + m]);
    }
#endif
  }
  for (unsigned int m = 0; m < iNetTop.NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor * neigh_proc_p =
        iNetTop.NeighbouringProcs[m];

    if (neigh_proc_p->Rank > iNetTop.LocalRank)
    {
#ifndef NOMPI
      err = MPI_Wait(&req[0][m], status);
#endif
    }
    else
    {
#ifndef NOMPI
      err = MPI_Wait(&req[0][iNetTop.NeighbouringProcs.size() + m], status);
#endif
      // Now we sort the situation so that each process has its own sites.
      for (int n = 0; n < neigh_proc_p->SharedFCount * 4; n += 4)
      {
        short int *f_data_p = &neigh_proc_p->SharedFLocation[n];

        short int l = f_data_p[3];
        f_data_p[0] += D3Q15::CX[l];
        f_data_p[1] += D3Q15::CY[l];
        f_data_p[2] += D3Q15::CZ[l];
        f_data_p[3] = D3Q15::INVERSEDIRECTIONS[l];
      }
    }
  }

  int f_count = bLocalLatDat.LocalFluidSites * D3Q15::NUMVECTORS;

  for (unsigned int m = 0; m < iNetTop.NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor *neigh_proc_p =
        iNetTop.NeighbouringProcs[m];

    for (int n = 0; n < neigh_proc_p->SharedFCount; n++)
    {
      // Get coordinates and direction of the distribution function to be sent to another process.
      short int *f_data_p = &neigh_proc_p->SharedFLocation[n * 4];
      short int i = f_data_p[0];
      short int j = f_data_p[1];
      short int k = f_data_p[2];
      short int l = f_data_p[3];

      // Get the fluid site number of site that will send data to another process.
      unsigned int site_map = *iGlobLatDat.GetSiteData(i, j, k);

      // Set f_id to the element in the send buffer that we put the updated
      // distribution functions in.
      bLocalLatDat.mFNeighbours[site_map * D3Q15::NUMVECTORS + l] = ++f_count;

      // Set the place where we put the received distribution functions, which is
      // f_new[number of fluid site that sends, inverse direction].
      neigh_proc_p->SharedFReceivingIndex[n] = site_map * D3Q15::NUMVECTORS
          + D3Q15::INVERSEDIRECTIONS[l];
    }
  }
  // neigh_prc->f_data was only set as a pointer to f_data, not allocated.  In this line, we
  // are freeing both of those.
  delete[] f_data;

  bm_time = hemelb::util::myClock() - seconds;
}

void Net::ReceiveFromNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat)
{
#ifndef NOMPI
  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor * neigh_proc_p =
        mNetworkTopology->NeighbouringProcs[m];

    err = MPI_Irecv(&bLocalLatDat.FOld[neigh_proc_p->FirstSharedF],
                    neigh_proc_p->SharedFCount, MPI_DOUBLE, neigh_proc_p->Rank,
                    10, MPI_COMM_WORLD, &req[0][m]);
  }
#endif
}

int Net::GetMachineCount()
{
  return mNetworkTopology->MachineCount;
}

void Net::SendToNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat)
{
#ifndef NOMPI
  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor *neigh_proc_p =
        mNetworkTopology->NeighbouringProcs[m];

    err = MPI_Isend(&bLocalLatDat.FNew[neigh_proc_p->FirstSharedF],
                    neigh_proc_p->SharedFCount, MPI_DOUBLE, neigh_proc_p->Rank,
                    10, MPI_COMM_WORLD,
                    &req[0][mNetworkTopology->NeighbouringProcs.size() + m]);
  }
#endif
}

void Net::UseDataFromNeighbouringProcs(hemelb::lb::LocalLatticeData &bLocalLatDat)
{
#ifndef NOMPI
  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    err = MPI_Wait(&req[0][m], status);
    err = MPI_Wait(&req[0][mNetworkTopology->NeighbouringProcs.size() + m],
                   status);
  }
#endif

  // Copy the distribution functions received from the neighbouring
  // processors into the destination buffer "f_new".
  for (int i = 0; i < mNetworkTopology->TotalSharedFs; i++)
  {
    bLocalLatDat.FNew[f_recv_iv[i]]
        = bLocalLatDat.FOld[mNetworkTopology->NeighbouringProcs[0]->FirstSharedF
            + i];
  }
}

Net::Net(hemelb::topology::NetworkTopology * iNetworkTopology,
         int &iArgumentCount,
         char* iArgumentList[])
{
#ifndef NOMPI
  int thread_level_provided;

  err = MPI_Init_thread(&iArgumentCount, &iArgumentList, MPI_THREAD_FUNNELED,
                        &thread_level_provided);
  err = MPI_Comm_size(MPI_COMM_WORLD, &iNetworkTopology->ProcessorCount);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &iNetworkTopology->LocalRank);

  if (iNetworkTopology->IsCurrentProcTheIOProc())
  {
    printf("thread_level_provided %i\n", thread_level_provided);
  }
#else
  iNetTop.ProcessorCount = 1;
  mRank = 0;
#endif

  mNetworkTopology = iNetworkTopology;
}

/*!
 Free the allocated data.
 */
Net::~Net()
{

#ifndef NOMPI
  err = MPI_Finalize();
#endif

  delete[] f_recv_iv;

#ifndef NOMPI
  for (int i = 0; i < COMMS_LEVELS; i++)
    delete[] req[i];
  delete[] req;
#endif

}

