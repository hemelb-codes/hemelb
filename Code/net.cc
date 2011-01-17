/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include "net.h"
#include "util/utilityFunctions.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

/*!
 This is called from the main function.  First function to deal with processors.
 The domain partitioning technique and the management of the
 buffers useful for the inter-processor communications are
 implemented in this function.  The domain decomposition is based
 on a graph growing partitioning technique.
 */
void Net::Initialise(hemelb::lb::GlobalLatticeData &iGlobLatDat,
                     hemelb::lb::LocalLatticeData* &bLocalLatDat)
{
  double seconds = hemelb::util::myClock();

  // Create a map between the two-level data representation and the 1D
  // compact one is created here.

  // This rank's site data.
  unsigned int
      *lThisRankSiteData =
          new unsigned int[mNetworkTopology->FluidSitesOnEachProcessor[mNetworkTopology->GetLocalRank()]];

  // Array of booleans to store whether any sites on a block are fluid
  // sites residing on this rank.
  bool *lBlockIsOnThisRank = new bool[iGlobLatDat.GetBlockCount()];
  // Initialise to false.
  for (int n = 0; n < iGlobLatDat.GetBlockCount(); n++)
  {
    lBlockIsOnThisRank[n] = false;
  }

  int lSiteIndexOnProc = 0;

  for (int lBlockNumber = 0; lBlockNumber < iGlobLatDat.GetBlockCount(); lBlockNumber++)
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
      if (mNetworkTopology->GetLocalRank()
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
  for (int n = 0; n < iGlobLatDat.GetBlockCount(); n++)
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

  mNetworkTopology->TotalSharedFs = 0; // shared SharedFCount within Net struct.

  lSiteIndexOnProc = 0;

  int n = -1;

  // Iterate over all blocks in site units
  for (int i = 0; i < iGlobLatDat.GetXSiteCount(); i
      += iGlobLatDat.GetBlockSize())
  {
    for (int j = 0; j < iGlobLatDat.GetYSiteCount(); j
        += iGlobLatDat.GetBlockSize())
    {
      for (int k = 0; k < iGlobLatDat.GetZSiteCount(); k
          += iGlobLatDat.GetBlockSize())
      {
        hemelb::lb::BlockData * map_block_p = &iGlobLatDat.Blocks[++n];

        if (map_block_p->site_data == NULL)
        {
          continue;
        }

        int m = -1;

        // Iterate over all sites within the current block.
        for (int site_i = i; site_i < i + iGlobLatDat.GetBlockSize(); site_i++)
        {
          for (int site_j = j; site_j < j + iGlobLatDat.GetBlockSize(); site_j++)
          {
            for (int site_k = k; site_k < k + iGlobLatDat.GetBlockSize(); site_k++)
            {
              m++;
              // If the site is not on this processor, continue.
              if (mNetworkTopology->GetLocalRank()
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
                // BIG_NUMBER2) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                // in lbmReadConfig in io.cc.
                if (proc_id_p == NULL || mNetworkTopology->GetLocalRank()
                    == (*proc_id_p) || *proc_id_p == (BIG_NUMBER2))
                {
                  continue;
                }

                lIsInnerSite = false;

                // The first time, net_neigh_procs = 0, so
                // the loop is not executed.
                bool flag = true;

                // Iterate over neighbouring processors until we find the one with the
                // neighbouring site on it.
                int lNeighbouringProcs =
                    mNetworkTopology->NeighbouringProcs.size();
                for (int mm = 0; mm < lNeighbouringProcs && flag; mm++)
                {
                  // Check whether the rank for a particular neighbour has already been
                  // used for this processor.  If it has, set flag to zero.
                  hemelb::topology::NeighbouringProcessor * neigh_proc_p =
                      mNetworkTopology->NeighbouringProcs[mm];

                  // If ProcessorRankForEachBlockSite is equal to a neigh_proc that has alredy been listed.
                  if (*proc_id_p == neigh_proc_p->Rank)
                  {
                    flag = false;
                    ++neigh_proc_p->SharedFCount;
                    ++mNetworkTopology->TotalSharedFs;
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
                  mNetworkTopology->NeighbouringProcs.push_back(lNewNeighbour);
                  ++mNetworkTopology->TotalSharedFs;
                }
              }

              // Set the collision type data. map_block site data is renumbered according to
              // fluid site numbers within a particular collision type.

              int l = -1;

              switch (iGlobLatDat.GetCollisionType(
                                                   lThisRankSiteData[lSiteIndexOnProc]))
              {
                case FLUID:
                  l = 0;
                  break;
                case EDGE:
                  l = 1;
                  break;
                case INLET:
                  l = 2;
                  break;
                case OUTLET:
                  l = 3;
                  break;
                case (INLET | EDGE):
                  l = 4;
                  break;
                case (OUTLET | EDGE):
                  l = 5;
                  break;
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

  bLocalLatDat
      = new hemelb::lb::LocalLatticeData(
                                         mNetworkTopology->FluidSitesOnEachProcessor[mNetworkTopology->GetLocalRank()],
                                         mNetworkTopology->TotalSharedFs);

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
  for (int n = 0; n < iGlobLatDat.GetBlockCount(); n++)
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

  short int *f_data = new short int[4 * mNetworkTopology->TotalSharedFs];

  // Allocate the index in which to put the distribution functions received from the other
  // process.
  f_recv_iv = new int[mNetworkTopology->TotalSharedFs];

  short int ** lSharedFLocationForEachProc =
      new short int*[mNetworkTopology->NeighbouringProcs.size()];

  // Reset to zero again.
  mNetworkTopology->TotalSharedFs = 0;

  // Set the remaining neighbouring processor data.
  for (unsigned int n = 0; n < mNetworkTopology->NeighbouringProcs.size(); n++)
  {
    // f_data compacted according to number of shared f_s on each process.
    // f_data will be set later.


    // Site co-ordinates for each of the shared distribution, and the number of
    // the corresponding direction vector. Array is 4 elements for each shared distribution.
    lSharedFLocationForEachProc[n] = &f_data[mNetworkTopology->TotalSharedFs
        << 2];

    // Pointing to a few things, but not setting any variables.
    // FirstSharedF points to start of shared_fs.
    mNetworkTopology->NeighbouringProcs[n]->FirstSharedF
        = bLocalLatDat->GetLocalFluidSiteCount() * D3Q15::NUMVECTORS + 1
            + mNetworkTopology->TotalSharedFs;

    mNetworkTopology->NeighbouringProcs[n]->SharedFReceivingIndex
        = &f_recv_iv[mNetworkTopology->TotalSharedFs];

    mNetworkTopology->TotalSharedFs
        += mNetworkTopology->NeighbouringProcs[n]->SharedFCount;
  }

  mNetworkTopology->NeighbourIndexFromProcRank
      = new short int[mNetworkTopology->GetProcessorCount()];

  for (int m = 0; m < mNetworkTopology->GetProcessorCount(); m++)
  {
    mNetworkTopology->NeighbourIndexFromProcRank[m] = -1;
  }
  // Get neigh_proc_index from ProcessorRankForEachBlockSite.
  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    mNetworkTopology->NeighbourIndexFromProcRank[mNetworkTopology->NeighbouringProcs[m]->Rank]
        = m;
  }

  InitialiseNeighbourLookup(bLocalLatDat, lSharedFLocationForEachProc,
                            lThisRankSiteData, iGlobLatDat);

  delete[] lThisRankSiteData;

  // point-to-point communications are performed to match data to be
  // sent to/receive from different partitions; in this way, the
  // communication of the locations of the interface-dependent fluid
  // sites and the identifiers of the distribution functions which
  // propagate to different partitions is avoided (only their values
  // will be communicated). It's here!

  // Allocate the request variable.
  req = new MPI_Request*[COMMS_LEVELS];

  for (int m = 0; m < COMMS_LEVELS; m++)
  {
    req[m] = new MPI_Request[2 * mNetworkTopology->GetProcessorCount()];
  }

  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor * neigh_proc_p =
        mNetworkTopology->NeighbouringProcs[m];

    // One way send receive.  The lower numbered mNetworkTopology->ProcessorCount send and the higher numbered ones receive.
    // It seems that, for each pair of processors, the lower numbered one ends up with its own
    // edge sites and directions stored and the higher numbered one ends up with those on the
    // other processor.
    if (neigh_proc_p->Rank > mNetworkTopology->GetLocalRank())
    {
      err = MPI_Isend(&lSharedFLocationForEachProc[m][0],
                      neigh_proc_p->SharedFCount * 4, MPI_SHORT,
                      neigh_proc_p->Rank, 10, MPI_COMM_WORLD, &req[0][m]);
    }
    else
    {
      err = MPI_Irecv(&lSharedFLocationForEachProc[m][0],
                      neigh_proc_p->SharedFCount * 4, MPI_SHORT,
                      neigh_proc_p->Rank, 10, MPI_COMM_WORLD,
                      &req[0][mNetworkTopology->NeighbouringProcs.size() + m]);
    }
  }
  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor * neigh_proc_p =
        mNetworkTopology->NeighbouringProcs[m];

    if (neigh_proc_p->Rank > mNetworkTopology->GetLocalRank())
    {
      err = MPI_Wait(&req[0][m], status);
    }
    else
    {
      err = MPI_Wait(&req[0][mNetworkTopology->NeighbouringProcs.size() + m],
                     status);

      // Now we sort the situation so that each process has its own sites.
      for (int n = 0; n < neigh_proc_p->SharedFCount * 4; n += 4)
      {
        short int *f_data_p = &lSharedFLocationForEachProc[m][n];

        short int l = f_data_p[3];
        f_data_p[0] += D3Q15::CX[l];
        f_data_p[1] += D3Q15::CY[l];
        f_data_p[2] += D3Q15::CZ[l];
        f_data_p[3] = D3Q15::INVERSEDIRECTIONS[l];
      }
    }
  }

  int f_count = bLocalLatDat->GetLocalFluidSiteCount() * D3Q15::NUMVECTORS;

  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    hemelb::topology::NeighbouringProcessor *neigh_proc_p =
        mNetworkTopology->NeighbouringProcs[m];

    for (int n = 0; n < neigh_proc_p->SharedFCount; n++)
    {
      // Get coordinates and direction of the distribution function to be sent to another process.
      short int *f_data_p = &lSharedFLocationForEachProc[m][n * 4];
      short int i = f_data_p[0];
      short int j = f_data_p[1];
      short int k = f_data_p[2];
      short int l = f_data_p[3];

      // Get the fluid site number of site that will send data to another process.
      unsigned int site_map = *iGlobLatDat.GetSiteData(i, j, k);

      // Set f_id to the element in the send buffer that we put the updated
      // distribution functions in.
      bLocalLatDat->SetNeighbourLocation(site_map, l, ++f_count);

      // Set the place where we put the received distribution functions, which is
      // f_new[number of fluid site that sends, inverse direction].
      neigh_proc_p->SharedFReceivingIndex[n] = site_map * D3Q15::NUMVECTORS
          + D3Q15::INVERSEDIRECTIONS[l];
    }
  }
  // neigh_prc->f_data was only set as a pointer to f_data, not allocated.  In this line, we
  // are freeing both of those.
  delete[] f_data;

  // Delete the array in which we kept the shared f locations. Don't delete subarrays - these
  // are pointers to elsewhere.
  delete[] lSharedFLocationForEachProc;

  bm_time = hemelb::util::myClock() - seconds;
}

void Net::InitialiseNeighbourLookup(hemelb::lb::LocalLatticeData* bLocalLatDat,
                                    short int ** bSharedFLocationForEachProc,
                                    const unsigned int * iSiteDataForThisRank,
                                    const hemelb::lb::GlobalLatticeData & iGlobLatDat)
{
  int n = -1;
  int lSiteIndexOnProc = 0;
  int * lFluidSitesHandledForEachProc =
      new int[mNetworkTopology->GetProcessorCount()];

  for (int ii = 0; ii < mNetworkTopology->GetProcessorCount(); ii++)
  {
    lFluidSitesHandledForEachProc[ii] = 0;
  }

  // Iterate over blocks in global co-ords.
  for (int i = 0; i < iGlobLatDat.GetXSiteCount(); i
      += iGlobLatDat.GetBlockSize())
  {
    for (int j = 0; j < iGlobLatDat.GetYSiteCount(); j
        += iGlobLatDat.GetBlockSize())
    {
      for (int k = 0; k < iGlobLatDat.GetZSiteCount(); k
          += iGlobLatDat.GetBlockSize())
      {
        n++;
        hemelb::lb::BlockData *map_block_p = &iGlobLatDat.Blocks[n];

        if (map_block_p->site_data == NULL)
        {
          continue;
        }

        int m = -1;

        // Iterate over sites within the block.
        for (int site_i = i; site_i < i + iGlobLatDat.GetBlockSize(); site_i++)
        {
          for (int site_j = j; site_j < j + iGlobLatDat.GetBlockSize(); site_j++)
          {
            for (int site_k = k; site_k < k + iGlobLatDat.GetBlockSize(); site_k++)
            {
              // If a site is not on this process, continue.
              m++;

              if (mNetworkTopology->GetLocalRank()
                  != map_block_p->ProcessorRankForEachBlockSite[m])
              {
                continue;
              }

              // Get site data, which is the number of the fluid site on this proc..
              unsigned int site_map = map_block_p->site_data[m];

              // Set neighbour location for the distribution component at the centre of
              // this site.
              bLocalLatDat->SetNeighbourLocation(site_map, 0, site_map
                  * D3Q15::NUMVECTORS + 0);

              for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
              {
                // Work out positions of neighbours.
                int neigh_i = site_i + D3Q15::CX[l];
                int neigh_j = site_j + D3Q15::CY[l];
                int neigh_k = site_k + D3Q15::CZ[l];

                // Get the id of the processor which the neighbouring site lies on.
                int *proc_id_p = iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i,
                                                                       neigh_j,
                                                                       neigh_k);

                if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
                {
                  // initialize f_id to the rubbish site.
                  bLocalLatDat->SetNeighbourLocation(
                                                     site_map,
                                                     l,
                                                     bLocalLatDat->GetLocalFluidSiteCount()
                                                         * D3Q15::NUMVECTORS);
                  continue;
                }
                // If on the same proc, set f_id of the
                // current site and direction to the
                // site and direction that it sends to.
                // If we check convergence, the data for
                // each site is split into that for the
                // current and previous cycles.
                else if (mNetworkTopology->GetLocalRank() == *proc_id_p)
                {

                  // Pointer to the neighbour.
                  const unsigned int *site_data_p =
                      iGlobLatDat.GetSiteData(neigh_i, neigh_j, neigh_k);

                  bLocalLatDat->SetNeighbourLocation(site_map, l, *site_data_p
                      * D3Q15::NUMVECTORS + l);

                  continue;
                }
                else
                {
                  short int neigh_proc_index =
                      mNetworkTopology->NeighbourIndexFromProcRank[*proc_id_p];

                  // This stores some coordinates.  We
                  // still need to know the site number.
                  // neigh_proc[ n ].f_data is now
                  // set as well, since this points to
                  // f_data.  Every process has data for
                  // its neighbours which say which sites
                  // on this process are shared with the
                  // neighbour.
                  short int
                      *f_data_p =
                          &bSharedFLocationForEachProc[neigh_proc_index][lFluidSitesHandledForEachProc[neigh_proc_index]
                              << 2];
                  f_data_p[0] = site_i;
                  f_data_p[1] = site_j;
                  f_data_p[2] = site_k;
                  f_data_p[3] = l;
                  ++lFluidSitesHandledForEachProc[neigh_proc_index];
                }
              }

              // This is used in Calculate BC in IO.
              bLocalLatDat->mSiteData[site_map]
                  = iSiteDataForThisRank[lSiteIndexOnProc];

              if (iGlobLatDat.GetCollisionType(
                                               bLocalLatDat->mSiteData[site_map])
                  & EDGE)
              {
                bLocalLatDat->SetWallNormal(
                                            site_map,
                                            iGlobLatDat.Blocks[n].wall_data[m].wall_nor);

                bLocalLatDat->SetDistanceToWall(
                                                site_map,
                                                iGlobLatDat.Blocks[n].wall_data[m].cut_dist);
              }
              else
              {
                double lBigDistance[3];
                for (unsigned int ii = 0; ii < 3; ii++)
                  lBigDistance[ii] = BIG_NUMBER;
                bLocalLatDat->SetWallNormal(site_map, lBigDistance);
              }
              ++lSiteIndexOnProc;
            }
          }
        }
      }
    }
  }

  delete[] lFluidSitesHandledForEachProc;
}

void Net::ReceiveFromNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat)
{
  int m = 0;

  for (std::vector<hemelb::topology::NeighbouringProcessor*>::iterator it =
      mNetworkTopology->NeighbouringProcs.begin(); it
      != mNetworkTopology->NeighbouringProcs.end(); ++it)
  {
    err = MPI_Irecv(&bLocalLatDat.FOld[ (*it)->FirstSharedF],
                     (*it)->SharedFCount, MPI_DOUBLE, (*it)->Rank, 10,
                    MPI_COMM_WORLD, &req[0][m]);
    ++m;
  }
}

void Net::SendToNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat)
{
  int m = 0;

  for (std::vector<hemelb::topology::NeighbouringProcessor*>::iterator it =
      mNetworkTopology->NeighbouringProcs.begin(); it
      != mNetworkTopology->NeighbouringProcs.end(); ++it)
  {
    err = MPI_Isend(&bLocalLatDat.FNew[ (*it)->FirstSharedF],
                     (*it)->SharedFCount, MPI_DOUBLE, (*it)->Rank, 10,
                    MPI_COMM_WORLD,
                    &req[0][mNetworkTopology->NeighbouringProcs.size() + m]);

    ++m;
  }
}

void Net::UseDataFromNeighbouringProcs(hemelb::lb::LocalLatticeData &bLocalLatDat)
{
  for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
  {
    err = MPI_Wait(&req[0][m], status);
    err = MPI_Wait(&req[0][mNetworkTopology->NeighbouringProcs.size() + m],
                   status);
  }

  // Copy the distribution functions received from the neighbouring
  // processors into the destination buffer "f_new".
  for (int i = 0; i < mNetworkTopology->TotalSharedFs; i++)
  {
    bLocalLatDat.FNew[f_recv_iv[i]]
        = bLocalLatDat.FOld[mNetworkTopology->NeighbouringProcs[0]->FirstSharedF
            + i];
  }
}

Net::Net(hemelb::topology::NetworkTopology * iNetworkTopology)
{
  mNetworkTopology = iNetworkTopology;
}

/*!
 Free the allocated data.
 */
Net::~Net()
{
  delete[] f_recv_iv;

  for (int i = 0; i < COMMS_LEVELS; i++)
    delete[] req[i];
  delete[] req;
}

