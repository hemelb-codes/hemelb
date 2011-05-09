/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "net/net.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace net
  {

    void Net::EnsureEnoughRequests(size_t count)
    {
      if (mRequests.size() < count)
      {
        size_t deficit = count - mRequests.size();
        for (unsigned int ii = 0; ii < deficit; ii++)
        {
          mRequests.push_back(MPI_Request());
          mStatuses.push_back(MPI_Status());
        }
      }
    }

    /*!
     This is called from the main function.  First function to deal with processors.
     The domain partitioning technique and the management of the
     buffers useful for the inter-processor communications are
     implemented in this function.  The domain decomposition is based
     on a graph growing partitioning technique.
     */
    site_t* Net::Initialise(geometry::LatticeData* bLatDat)
    {
      // Create a map between the two-level data representation and the 1D
      // compact one is created here.

      // This rank's site data.
      unsigned int
          *lThisRankSiteData =
              new unsigned int[mNetworkTopology->FluidSitesOnEachProcessor[mNetworkTopology->GetLocalRank()]];

      GetThisRankSiteData(bLatDat, lThisRankSiteData);

      // The numbers of inter- and intra-machine neighbouring processors,
      // interface-dependent and independent fluid sites and shared
      // distribution functions of the reference processor are calculated
      // here.  neigh_proc is a static array that is declared in config.h.

      CountCollisionTypes(bLatDat, lThisRankSiteData);

      // the precise interface-dependent data (interface-dependent fluid
      // site locations and identifiers of the distribution functions
      // streamed between different partitions) are collected and the
      // buffers needed for the communications are set from here

      site_t* f_data = new site_t[4 * mNetworkTopology->TotalSharedFs];

      // Allocate the index in which to put the distribution functions received from the other
      // process.
      site_t* f_recv_iv = new site_t[mNetworkTopology->TotalSharedFs];

      site_t** lSharedFLocationForEachProc =
          new site_t*[mNetworkTopology->NeighbouringProcs.size()];

      // Reset to zero again.
      mNetworkTopology->TotalSharedFs = 0;

      // Set the remaining neighbouring processor data.
      for (size_t n = 0; n < mNetworkTopology->NeighbouringProcs.size(); n++)
      {
        // f_data compacted according to number of shared f_s on each process.
        // f_data will be set later.


        // Site co-ordinates for each of the shared distribution, and the number of
        // the corresponding direction vector. Array is 4 elements for each shared distribution.
        lSharedFLocationForEachProc[n] = &f_data[mNetworkTopology->TotalSharedFs << 2];

        // Pointing to a few things, but not setting any variables.
        // FirstSharedF points to start of shared_fs.
        mNetworkTopology->NeighbouringProcs[n].FirstSharedF = bLatDat->GetLocalFluidSiteCount()
            * D3Q15::NUMVECTORS + 1 + mNetworkTopology->TotalSharedFs;

        mNetworkTopology->TotalSharedFs += mNetworkTopology->NeighbouringProcs[n].SharedFCount;
      }

      mNetworkTopology->NeighbourIndexFromProcRank
          = new proc_t[mNetworkTopology->GetProcessorCount()];

      for (proc_t m = 0; m < mNetworkTopology->GetProcessorCount(); m++)
      {
        mNetworkTopology->NeighbourIndexFromProcRank[m] = -1;
      }
      // Get neigh_proc_index from ProcessorRankForEachBlockSite.
      for (proc_t m = 0; m < (proc_t) mNetworkTopology->NeighbouringProcs.size(); m++)
      {
        mNetworkTopology->NeighbourIndexFromProcRank[mNetworkTopology->NeighbouringProcs[m].Rank]
            = m;
      }

      {
        site_t** SharedLocationPerProcByNeighbourId =
            new site_t*[mNetworkTopology->GetProcessorCount()];

        for (proc_t ii = 0; ii < mNetworkTopology->GetProcessorCount(); ++ii)
        {
          if (mNetworkTopology->NeighbourIndexFromProcRank[ii] >= 0)
          {
            SharedLocationPerProcByNeighbourId[ii]
                = lSharedFLocationForEachProc[mNetworkTopology->NeighbourIndexFromProcRank[ii]];
          }
        }

        bLatDat->InitialiseNeighbourLookup(SharedLocationPerProcByNeighbourId,
                                           mNetworkTopology->GetLocalRank(),
                                           lThisRankSiteData);

        delete[] SharedLocationPerProcByNeighbourId;
      }

      delete[] lThisRankSiteData;

      InitialisePointToPointComms(lSharedFLocationForEachProc);

      site_t f_count = bLatDat->GetLocalFluidSiteCount() * D3Q15::NUMVECTORS;

      site_t sharedSitesSeen = 0;

      for (size_t m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
      {
        hemelb::topology::NeighbouringProcessor *neigh_proc_p =
            &mNetworkTopology->NeighbouringProcs[m];

        for (site_t n = 0; n < neigh_proc_p->SharedFCount; n++)
        {
          // Get coordinates and direction of the distribution function to be sent to another process.
          site_t* f_data_p = &lSharedFLocationForEachProc[m][n * 4];
          site_t i = f_data_p[0];
          site_t j = f_data_p[1];
          site_t k = f_data_p[2];
          site_t l = f_data_p[3];

          // Correct so that each process has the correct coordinates.
          if (neigh_proc_p->Rank < mNetworkTopology->GetLocalRank())
          {
            i += D3Q15::CX[l];
            j += D3Q15::CY[l];
            k += D3Q15::CZ[l];
            l = D3Q15::INVERSEDIRECTIONS[l];
          }

          // Get the fluid site number of site that will send data to another process.
          site_t contigSiteId = bLatDat->GetContiguousSiteId(i, j, k);

          // Set f_id to the element in the send buffer that we put the updated
          // distribution functions in.
          bLatDat->SetNeighbourLocation(contigSiteId, (unsigned int) l, ++f_count);

          // Set the place where we put the received distribution functions, which is
          // f_new[number of fluid site that sends, inverse direction].
          f_recv_iv[sharedSitesSeen] = contigSiteId * D3Q15::NUMVECTORS
              + D3Q15::INVERSEDIRECTIONS[l];
          ++sharedSitesSeen;
        }
      }
      // neigh_prc->f_data was only set as a pointer to f_data, not allocated.  In this line, we
      // are freeing both of those.
      delete[] f_data;

      // Delete the array in which we kept the shared f locations. Don't delete subarrays - these
      // are pointers to elsewhere.
      delete[] lSharedFLocationForEachProc;

      return f_recv_iv;
    }

    void Net::InitialisePointToPointComms(site_t** &lSharedFLocationForEachProc)
    {
      EnsureEnoughRequests(mNetworkTopology->NeighbouringProcs.size());

      // point-to-point communications are performed to match data to be
      // sent to/receive from different partitions; in this way, the
      // communication of the locations of the interface-dependent fluid
      // sites and the identifiers of the distribution functions which
      // propagate to different partitions is avoided (only their values
      // will be communicated). It's here!
      // Allocate the request variable.
      for (size_t m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
      {
        hemelb::topology::NeighbouringProcessor *neigh_proc_p =
            &mNetworkTopology->NeighbouringProcs[m];
        // One way send receive.  The lower numbered mNetworkTopology->ProcessorCount send and the higher numbered ones receive.
        // It seems that, for each pair of processors, the lower numbered one ends up with its own
        // edge sites and directions stored and the higher numbered one ends up with those on the
        // other processor.
        if (neigh_proc_p->Rank > mNetworkTopology->GetLocalRank())
        {
          MPI_Isend(&lSharedFLocationForEachProc[m][0],
                    (int) neigh_proc_p->SharedFCount * 4,
                    MpiDataType<site_t>(),
                    neigh_proc_p->Rank,
                    10,
                    MPI_COMM_WORLD,
                    &mRequests[m]);
        }
        else
        {
          MPI_Irecv(&lSharedFLocationForEachProc[m][0],
                    (int) neigh_proc_p->SharedFCount * 4,
                    MpiDataType<site_t>(),
                    neigh_proc_p->Rank,
                    10,
                    MPI_COMM_WORLD,
                    &mRequests[m]);
        }
      }

      MPI_Waitall((int) mNetworkTopology->NeighbouringProcs.size(), &mRequests[0], &mStatuses[0]);
    }

    void Net::GetThisRankSiteData(const geometry::LatticeData* iLatDat,
                                  unsigned int* &bThisRankSiteData)
    {
      // Array of booleans to store whether any sites on a block are fluid
      // sites residing on this rank.
      bool *lBlockIsOnThisRank = new bool[iLatDat->GetBlockCount()];
      // Initialise to false.
      for (site_t n = 0; n < iLatDat->GetBlockCount(); n++)
      {
        lBlockIsOnThisRank[n] = false;
      }

      int lSiteIndexOnProc = 0;

      for (site_t lBlockNumber = 0; lBlockNumber < iLatDat->GetBlockCount(); lBlockNumber++)
      {
        geometry::LatticeData::BlockData* lCurrentDataBlock = iLatDat->GetBlock(lBlockNumber);

        // If we are in a block of solids, move to the next block.
        if (lCurrentDataBlock->site_data == NULL)
        {
          continue;
        }

        // If we have some fluid sites, point to mProcessorsForEachBlock and map_block.
        geometry::LatticeData::BlockData * proc_block_p = iLatDat->GetBlock(lBlockNumber);

        // lCurrentDataBlock.site_data is set to the fluid site identifier on this rank or (1U << 31U) if a site is solid
        // or not on this rank.  site_data is indexed by fluid site identifier and set to the site_data.
        for (site_t lSiteIndexWithinBlock = 0; lSiteIndexWithinBlock
            < iLatDat->GetSitesPerBlockVolumeUnit(); lSiteIndexWithinBlock++)
        {
          if (mNetworkTopology->GetLocalRank()
              == proc_block_p->ProcessorRankForEachBlockSite[lSiteIndexWithinBlock])
          {
            // If the current site is non-solid, copy the site data into the array for
            // this rank (in the whole-processor location), then set the site data
            // for this site within the current block to be the site index over the whole
            // processor.
            if ( (lCurrentDataBlock->site_data[lSiteIndexWithinBlock] & SITE_TYPE_MASK)
                != geometry::LatticeData::SOLID_TYPE)
            {
              bThisRankSiteData[lSiteIndexOnProc]
                  = lCurrentDataBlock->site_data[lSiteIndexWithinBlock];
              lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = lSiteIndexOnProc;
              ++lSiteIndexOnProc;
            }
            else
            {
              // If this is a solid, set the site data on the current block to
              // some massive value.
              lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = BIG_NUMBER3;
            }
            // Set the array to notify that the current block has sites on this
            // rank.
            lBlockIsOnThisRank[lBlockNumber] = true;
          }
          // If this site is not on the current processor, set its whole processor
          // index within the per-block store to a nonsense value.
          else
          {
            lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = BIG_NUMBER3;
          }
        }
      }

      // If we are in a block of solids, we set map_block[n].site_data to NULL.
      for (site_t n = 0; n < iLatDat->GetBlockCount(); n++)
      {
        if (lBlockIsOnThisRank[n])
        {
          continue;
        }

        if (iLatDat->GetBlock(n)->site_data != NULL)
        {
          delete[] iLatDat->GetBlock(n)->site_data;
          iLatDat->GetBlock(n)->site_data = NULL;
        }

        if (iLatDat->GetBlock(n)->wall_data != NULL)
        {
          delete[] iLatDat->GetBlock(n)->wall_data;
          iLatDat->GetBlock(n)->wall_data = NULL;
        }
      }
      delete[] lBlockIsOnThisRank;
    }

    void Net::CountCollisionTypes(geometry::LatticeData* bLatDat,
                                  const unsigned int * lThisRankSiteData)
    {
      site_t innerSites = 0;

      site_t interCollisions[COLLISION_TYPES];
      site_t innerCollisions[COLLISION_TYPES];

      for (unsigned int m = 0; m < COLLISION_TYPES; m++)
      {
        interCollisions[m] = 0;
        innerCollisions[m] = 0;
      }

      mNetworkTopology->TotalSharedFs = 0;

      int lSiteIndexOnProc = 0;

      site_t n = -1;

      // Iterate over all blocks in site units
      for (site_t i = 0; i < bLatDat->GetXSiteCount(); i += bLatDat->GetBlockSize())
      {
        for (site_t j = 0; j < bLatDat->GetYSiteCount(); j += bLatDat->GetBlockSize())
        {
          for (site_t k = 0; k < bLatDat->GetZSiteCount(); k += bLatDat->GetBlockSize())
          {
            geometry::LatticeData::BlockData * map_block_p = bLatDat->GetBlock(++n);

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            site_t m = -1;

            // Iterate over all sites within the current block.
            for (site_t site_i = i; site_i < i + bLatDat->GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + bLatDat->GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + bLatDat->GetBlockSize(); site_k++)
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
                    site_t neigh_i = site_i + D3Q15::CX[l];
                    site_t neigh_j = site_j + D3Q15::CY[l];
                    site_t neigh_k = site_k + D3Q15::CZ[l];

                    if (!bLatDat->IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      continue;
                    }

                    // Find the processor Id for that neighbour.
                    const proc_t* proc_id_p = bLatDat->GetProcIdFromGlobalCoords(neigh_i,
                                                                                 neigh_j,
                                                                                 neigh_k);

                    // Move on if the neighbour is in a block of solids (in which case
                    // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                    // BIG_NUMBER2) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                    // in lbmReadConfig in io.cc.
                    if (proc_id_p == NULL || mNetworkTopology->GetLocalRank() == *proc_id_p
                        || *proc_id_p == (BIG_NUMBER2))
                    {
                      continue;
                    }

                    lIsInnerSite = false;

                    // The first time, net_neigh_procs = 0, so
                    // the loop is not executed.
                    bool flag = true;

                    // Iterate over neighbouring processors until we find the one with the
                    // neighbouring site on it.
                    proc_t lNeighbouringProcs = (proc_t) mNetworkTopology->NeighbouringProcs.size();
                    for (proc_t mm = 0; mm < lNeighbouringProcs && flag; mm++)
                    {
                      // Check whether the rank for a particular neighbour has already been
                      // used for this processor.  If it has, set flag to zero.
                      hemelb::topology::NeighbouringProcessor* neigh_proc_p =
                          &mNetworkTopology->NeighbouringProcs[mm];

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
                      hemelb::topology::NeighbouringProcessor lNewNeighbour;
                      lNewNeighbour.SharedFCount = 1;
                      lNewNeighbour.Rank = *proc_id_p;
                      mNetworkTopology->NeighbouringProcs.push_back(lNewNeighbour);
                      ++mNetworkTopology->TotalSharedFs;
                    }
                  }

                  // Set the collision type data. map_block site data is renumbered according to
                  // fluid site numbers within a particular collision type.

                  int l = -1;

                  switch (bLatDat->GetCollisionType(lThisRankSiteData[lSiteIndexOnProc]))
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
                    ++innerSites;

                    if (l == 0)
                    {
                      map_block_p->site_data[m] = (unsigned int) innerCollisions[l];
                    }

                    ++innerCollisions[l];
                  }
                  else
                  {
                    ++interCollisions[l];
                  }
                }
              }
            }
          }
        }
      }

      site_t collision_offset[2][COLLISION_TYPES];
      // Calculate the number of each type of collision.
      collision_offset[0][0] = 0;

      for (unsigned int l = 1; l < COLLISION_TYPES; l++)
      {
        collision_offset[0][l] = collision_offset[0][l - 1] + innerCollisions[l - 1];
      }
      collision_offset[1][0] = innerSites;
      for (unsigned int l = 1; l < COLLISION_TYPES; l++)
      {
        collision_offset[1][l] = collision_offset[1][l - 1] + interCollisions[l - 1];
      }

      lSiteIndexOnProc = 0;

      site_t innerColsPassed[COLLISION_TYPES];
      site_t interColsPassed[COLLISION_TYPES];

      for (unsigned int ii = 0; ii < COLLISION_TYPES; ++ii)
      {
        innerColsPassed[ii] = 0;
        interColsPassed[ii] = 0;
      }

      // Iterate over blocks
      n = -1;

      // Iterate over all blocks in site units
      for (site_t i = 0; i < bLatDat->GetXSiteCount(); i += bLatDat->GetBlockSize())
      {
        for (site_t j = 0; j < bLatDat->GetYSiteCount(); j += bLatDat->GetBlockSize())
        {
          for (site_t k = 0; k < bLatDat->GetZSiteCount(); k += bLatDat->GetBlockSize())
          {
            geometry::LatticeData::BlockData * map_block_p = bLatDat->GetBlock(++n);

            // If we are in a block of solids, continue.
            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            // Iterate over sites within the block.
            site_t m = -1;

            // Iterate over all sites within the current block.
            for (site_t site_i = i; site_i < i + bLatDat->GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + bLatDat->GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + bLatDat->GetBlockSize(); site_k++)
                {
                  m++;

                  unsigned int *site_data_p = &map_block_p->site_data[m];

                  // If the site is solid, continue.
                  if (*site_data_p & BIG_NUMBER3)
                  {
                    continue;
                  }

                  int l = -1;

                  switch (bLatDat->GetCollisionType(lThisRankSiteData[lSiteIndexOnProc]))
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

                  bool lIsInnerSite = true;

                  // Iterate over all direction vectors.
                  for (unsigned int q = 1; q < D3Q15::NUMVECTORS; q++)
                  {
                    // Find the neighbour site co-ords in this direction.
                    site_t neigh_i = site_i + D3Q15::CX[q];
                    site_t neigh_j = site_j + D3Q15::CY[q];
                    site_t neigh_k = site_k + D3Q15::CZ[q];

                    if (!bLatDat->IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      continue;
                    }

                    // Find the processor Id for that neighbour.
                    const proc_t* proc_id_p = bLatDat->GetProcIdFromGlobalCoords(neigh_i,
                                                                                 neigh_j,
                                                                                 neigh_k);

                    // Move on if the neighbour is in a block of solids (in which case
                    // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                    // BIG_NUMBER2) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                    // in lbmReadConfig in io.cc.
                    if (proc_id_p == NULL || (int) mNetworkTopology->GetLocalRank() == *proc_id_p
                        || *proc_id_p == (BIG_NUMBER2))
                    {
                      continue;
                    }

                    lIsInnerSite = false;
                    break;
                  }

                  if (lIsInnerSite)
                  {
                    if (l != 0)
                    {
                      *site_data_p = (unsigned int) (collision_offset[0][l] + innerColsPassed[l]);
                      ++innerColsPassed[l];
                    }
                  }
                  else
                  {
                    *site_data_p = (unsigned int) (collision_offset[1][l] + interColsPassed[l]);
                    ++interColsPassed[l];
                  }
                }
              }
            }
          }
        }
      }

      bLatDat->SetSiteCounts(innerSites,
                             interCollisions,
                             innerCollisions,
                             mNetworkTopology->TotalSharedFs);
    }

    void Net::Receive()
    {
      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();

      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); ++it)
      {
        MPI_Irecv(it->second.PointerList.front(),
                  1,
                  it->second.Type,
                  it->first,
                  10,
                  MPI_COMM_WORLD,
                  &mRequests[m]);
        ++m;
      }
    }

    void Net::Send()
    {
      // Make sure the datatypes have been created.
      EnsurePreparedToSendReceive();

      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); ++it)
      {
        MPI_Isend(it->second.PointerList.front(),
                  1,
                  it->second.Type,
                  it->first,
                  10,
                  MPI_COMM_WORLD,
                  &mRequests[mReceiveProcessorComms.size() + m]);

        ++m;
      }
    }

    void Net::Wait()
    {
      MPI_Waitall((int) (mSendProcessorComms.size() + mReceiveProcessorComms.size()),
                  &mRequests[0],
                  &mStatuses[0]);

      sendReceivePrepped = false;

      for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); it++)
      {
        MPI_Type_free(&it->second.Type);
      }
      mReceiveProcessorComms.clear();

      for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); it++)
      {
        MPI_Type_free(&it->second.Type);
      }
      mSendProcessorComms.clear();
    }

    // Helper function to get the ProcessorCommunications object, and create it if it doesn't exist yet.
    Net::ProcComms* Net::GetProcComms(proc_t iRank, bool iIsSend)
    {
      std::map<proc_t, ProcComms>* lMap = iIsSend
        ? &mSendProcessorComms
        : &mReceiveProcessorComms;

      std::map<proc_t, ProcComms>::iterator lValue = lMap->find(iRank);

      if (lValue == lMap->end())
      {
        ProcComms lRet;
        lMap->insert(std::pair<proc_t, ProcComms>(iRank, lRet));
        return &lMap->find(iRank)->second;
      }
      else
      {
        return & (lValue ->second);
      }
    }

//    // Helper functions to add ints to the list.
//    void Net::AddToList(int* iNew, int iLength, ProcComms *bMetaData)
//    {
//      bMetaData->PointerList.push_back(iNew);
//      bMetaData->LengthList.push_back(iLength);
//      bMetaData->TypeList.push_back(MpiDataType(iNew));
//    }
//
//    // Helper functions to add doubles to the list.
//    void Net::AddToList(double* iNew, int iLength, ProcComms *bMetaData)
//    {
//      bMetaData->PointerList.push_back(iNew);
//      bMetaData->LengthList.push_back(iLength);
//      bMetaData->TypeList.push_back(MpiDataType(iNew));
//    }
//
//    // Helper functions to add floats to the list.
//    void Net::AddToList(float* iNew, int iLength, ProcComms *bMetaData)
//    {
//      bMetaData->PointerList.push_back(iNew);
//      bMetaData->LengthList.push_back(iLength);
//      bMetaData->TypeList.push_back(MpiDataType(iNew));
//    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void Net::EnsurePreparedToSendReceive()
    {
      if (sendReceivePrepped)
      {
        return;
      }

      for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); it++)
      {
        CreateMPIType(& (*it).second);
      }

      for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); it++)
      {
        CreateMPIType(& (*it).second);
      }

      EnsureEnoughRequests(mReceiveProcessorComms.size() + mSendProcessorComms.size());

      sendReceivePrepped = true;
    }

    // Helper function to create a MPI derived datatype given a list of pointers, types and lengths.
    void Net::CreateMPIType(ProcComms *iMetaData)
    {
      MPI_Aint* displacements = new MPI_Aint[iMetaData->PointerList.size()];
      int* lengths = new int[iMetaData->PointerList.size()];
      MPI_Datatype* types = new MPI_Datatype[iMetaData->PointerList.size()];

      int lLocation = 0;

      for (std::vector<void*>::const_iterator it = iMetaData->PointerList.begin(); it
          != iMetaData->PointerList.end(); it++)
      {
        MPI_Get_address(*it, &displacements[lLocation]);
        ++lLocation;
      }

      for (int ii = (int) iMetaData->PointerList.size() - 1; ii >= 0; ii--)
      {
        displacements[ii] -= displacements[0];

        lengths[ii] = iMetaData->LengthList[ii];
        types[ii] = iMetaData->TypeList[ii];
      }

      // Create the type and commit it.
      MPI_Type_create_struct((int) iMetaData->PointerList.size(),
                             lengths,
                             displacements,
                             types,
                             &iMetaData->Type);
      MPI_Type_commit(&iMetaData->Type);

      delete[] displacements;
      delete[] lengths;
      delete[] types;
    }

    Net::Net(hemelb::topology::NetworkTopology * iNetworkTopology)
    {
      mNetworkTopology = iNetworkTopology;
      sendReceivePrepped = false;
    }

    /*!
     Free the allocated data.
     */
    Net::~Net()
    {
      if (sendReceivePrepped)
      {
        for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
            != mSendProcessorComms.end(); it++)
        {
          MPI_Type_free(&it->second.Type);
        }

        for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
            != mReceiveProcessorComms.end(); it++)
        {
          MPI_Type_free(&it->second.Type);
        }

      }
    }

  }
}
