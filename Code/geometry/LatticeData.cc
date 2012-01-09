#include <map>
#include <limits>

#include "debug/Debugger.h"
#include "log/Logger.h"
#include "topology/NetworkTopology.h"
#include "geometry/LatticeData.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    LatticeData::LatticeData()
    {
    }

    LatticeData::~LatticeData()
    {
    }

    LatticeData* LatticeData::Load(const bool reserveSteeringCore,
                                   std::string& dataFilePath,
                                   reporting::Timers &timings)
    {
      // Use a reader to read in the file.
      hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Loading file and decomposing geometry.");

      GeometryReadResult readGeometryData;
      GeometryReader reader(reserveSteeringCore, readGeometryData);
      reader.LoadAndDecompose(dataFilePath, timings);

      // Create a new lattice based on that info and return it.
      return new LatticeData(readGeometryData);
    }

    LatticeData::LatticeData(const GeometryReadResult& readResult)
    {
      SetBasicDetails(readResult.blocks,
                      readResult.blockSize,
                      readResult.voxelSize,
                      readResult.origin);

      ProcessReadSites(readResult);
      CollectFluidSiteDistribution();
      CollectGlobalSiteExtrema();

      InitialiseNeighbourStuff();

      CleanEmptyBlocks();
    }

    void LatticeData::SetBasicDetails(util::Vector3D<site_t> blocksIn,
                                      site_t blockSizeIn,
                                      distribn_t voxelSizeIn,
                                      util::Vector3D<distribn_t> originIn)
    {
      blockCounts = blocksIn;
      blockSize = blockSizeIn;
      voxelSize = voxelSizeIn;
      origin = originIn;
      sites = blocksIn * blockSize;
      sitesPerBlockVolumeUnit = blockSize * blockSize * blockSize;
      // A shift value we'll need later = log_2(block_size)
      site_t i = blockSize;
      log2BlockSize = 0;
      while (i > 1)
      {
        i >>= 1;
        ++log2BlockSize;
      }
      blockCount = blockCounts.x * blockCounts.y * blockCounts.z;
    }

    void LatticeData::ProcessReadSites(const GeometryReadResult & readResult)
    {
      Blocks.resize(GetBlockCount());

      totalSharedFs = 0;

      std::vector<SiteData> interSiteData[COLLISION_TYPES];
      std::vector<SiteData> intraSiteData[COLLISION_TYPES];
      std::vector<site_t> interBlockNumber[COLLISION_TYPES];
      std::vector<site_t> intraBlockNumber[COLLISION_TYPES];
      std::vector<site_t> interSiteNumber[COLLISION_TYPES];
      std::vector<site_t> intraSiteNumber[COLLISION_TYPES];
      std::vector<util::Vector3D<double> > interWallNormals[COLLISION_TYPES];
      std::vector<util::Vector3D<double> > intraWallNormals[COLLISION_TYPES];
      std::vector<double> interWallDistance[COLLISION_TYPES];
      std::vector<double> intraWallDistance[COLLISION_TYPES];

      proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();
      // Iterate over all blocks in site units
      for (site_t blockI = 0; blockI < blockCounts.x; blockI++)
      {
        for (site_t blockJ = 0; blockJ < blockCounts.y; blockJ++)
        {
          for (site_t blockK = 0; blockK < blockCounts.z; blockK++)
          {
            site_t blockId = GetBlockIdFromBlockCoords(blockI, blockJ, blockK);
            const BlockReadResult & blockReadIn = readResult.Blocks[blockId];
            if (blockReadIn.Sites.size() == 0)
            {
              continue;
            }
            // Iterate over all sites within the current block.
            for (site_t localSiteI = 0; localSiteI < readResult.blockSize; localSiteI++)
            {
              for (site_t localSiteJ = 0; localSiteJ < readResult.blockSize; localSiteJ++)
              {
                for (site_t localSiteK = 0; localSiteK < readResult.blockSize; localSiteK++)
                {
                  site_t localSiteId = readResult.GetSiteIdFromSiteCoordinates(localSiteI,
                                                                               localSiteJ,
                                                                               localSiteK);

                  if (Blocks[blockId].localContiguousIndex.size() == 0)
                  {
                    Blocks[blockId] = BlockData(GetSitesPerBlockVolumeUnit());
                  }

                  Blocks[blockId].processorRankForEachBlockSite[localSiteId]
                      = blockReadIn.Sites[localSiteId].targetProcessor;

                  // If the site is not on this processor, continue.
                  if (localRank != blockReadIn.Sites[localSiteId].targetProcessor)
                  {
                    continue;
                  }
                  bool lIsInnerSite = true;
                  // Iterate over all direction vectors.
                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Find the neighbour site co-ords in this direction.
                    site_t neighGlobalI = blockI * readResult.blockSize + localSiteI + D3Q15::CX[l];
                    site_t neighGlobalJ = blockJ * readResult.blockSize + localSiteJ + D3Q15::CY[l];
                    site_t neighGlobalK = blockK * readResult.blockSize + localSiteK + D3Q15::CZ[l];
                    if (neighGlobalI < 0 || neighGlobalJ < 0 || neighGlobalK < 0 || neighGlobalI
                        >= readResult.blocks.x * readResult.blockSize || neighGlobalJ
                        >= readResult.blocks.y * readResult.blockSize || neighGlobalK
                        >= readResult.blocks.z * readResult.blockSize)
                    {
                      continue;
                    }
                    // ... (that is actually being simulated and not a solid)...
                    site_t neighbourBlockI = neighGlobalI / readResult.blockSize;
                    site_t neighbourBlockJ = neighGlobalJ / readResult.blockSize;
                    site_t neighbourBlockK = neighGlobalK / readResult.blockSize;
                    site_t neighbourSiteI = neighGlobalI % readResult.blockSize;
                    site_t neighbourSiteJ = neighGlobalJ % readResult.blockSize;
                    site_t neighbourSiteK = neighGlobalK % readResult.blockSize;
                    site_t neighbourBlockId =
                        readResult.GetBlockIdFromBlockCoordinates(neighbourBlockI,
                                                                  neighbourBlockJ,
                                                                  neighbourBlockK);
                    // Move on if the neighbour is in a block of solids (in which case
                    // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                    // BIG_NUMBER2) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                    // in lbmReadConfig in io.cc.
                    if (readResult.Blocks[neighbourBlockId].Sites.size() == 0)
                    {
                      continue;
                    }
                    site_t neighbourSiteId =
                        readResult.GetSiteIdFromSiteCoordinates(neighbourSiteI,
                                                                neighbourSiteJ,
                                                                neighbourSiteK);
                    proc_t neighbourProc =
                        readResult.Blocks[neighbourBlockId].Sites[neighbourSiteId].targetProcessor;
                    if (neighbourProc == BIG_NUMBER2 || localRank == neighbourProc)
                    {
                      continue;
                    }
                    lIsInnerSite = false;
                    totalSharedFs++;

                    // The first time, net_neigh_procs = 0, so
                    // the loop is not executed.
                    bool flag = true;
                    // Iterate over neighbouring processors until we find the one with the
                    // neighbouring site on it.
                    proc_t lNeighbouringProcs = (proc_t) (neighbouringProcs.size());
                    for (proc_t mm = 0; mm < lNeighbouringProcs && flag; mm++)
                    {
                      // Check whether the rank for a particular neighbour has already been
                      // used for this processor.  If it has, set flag to zero.
                      NeighbouringProcessor *neigh_proc_p = &neighbouringProcs[mm];
                      // If ProcessorRankForEachBlockSite is equal to a neigh_proc that has alredy been listed.
                      if (neighbourProc == neigh_proc_p->Rank)
                      {
                        flag = false;
                        ++neigh_proc_p->SharedFCount;
                        break;
                      }
                    }

                    // If flag is 1, we didn't find a neighbour-proc with the neighbour-site on it
                    // so we need a new neighbouring processor.
                    if (flag)
                    {
                      // Store rank of neighbour in >neigh_proc[neigh_procs]
                      NeighbouringProcessor lNewNeighbour;
                      lNewNeighbour.SharedFCount = 1;
                      lNewNeighbour.Rank = neighbourProc;
                      neighbouringProcs.push_back(lNewNeighbour);
                    }
                  }

                  // Set the collision type data. map_block site data is renumbered according to
                  // fluid site numbers within a particular collision type.
                  int l = -1;
                  switch (blockReadIn.Sites[localSiteId].siteData.GetCollisionType())
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

                  if (lIsInnerSite)
                  {
                    intraBlockNumber[l].push_back(blockId);
                    intraSiteNumber[l].push_back(localSiteId);
                    intraSiteData[l].push_back(blockReadIn.Sites[localSiteId].siteData);
                    intraWallNormals[l].push_back(blockReadIn.Sites[localSiteId].wallNormal);
                    for (Direction direction = 1; direction < D3Q15::NUMVECTORS; direction++)
                    {
                      intraWallDistance[l].push_back(blockReadIn.Sites[localSiteId].cutDistance[direction
                          - 1]);
                    }
                  }
                  else
                  {
                    interBlockNumber[l].push_back(blockId);
                    interSiteNumber[l].push_back(localSiteId);
                    interSiteData[l].push_back(blockReadIn.Sites[localSiteId].siteData);
                    interWallNormals[l].push_back(blockReadIn.Sites[localSiteId].wallNormal);
                    for (Direction direction = 1; direction < D3Q15::NUMVECTORS; direction++)
                    {
                      interWallDistance[l].push_back(blockReadIn.Sites[localSiteId].cutDistance[direction
                          - 1]);
                    }
                  }
                }
              }
            }
          }
        }
      }

      PopulateWithReadData(intraBlockNumber,
                           intraSiteNumber,
                           intraSiteData,
                           intraWallNormals,
                           intraWallDistance,
                           interBlockNumber,
                           interSiteNumber,
                           interSiteData,
                           interWallNormals,
                           interWallDistance);
    }

    void LatticeData::PopulateWithReadData(const std::vector<site_t> intraBlockNumbers[COLLISION_TYPES],
                                           const std::vector<site_t> intraSiteNumbers[COLLISION_TYPES],
                                           const std::vector<SiteData> intraSiteData[COLLISION_TYPES],
                                           const std::vector<util::Vector3D<double> > intraWallNormals[COLLISION_TYPES],
                                           const std::vector<double> intraWallDistance[COLLISION_TYPES],
                                           const std::vector<site_t> interBlockNumbers[COLLISION_TYPES],
                                           const std::vector<site_t> interSiteNumbers[COLLISION_TYPES],
                                           const std::vector<SiteData> interSiteData[COLLISION_TYPES],
                                           const std::vector<util::Vector3D<double> > interWallNormals[COLLISION_TYPES],
                                           const std::vector<double> interWallDistance[COLLISION_TYPES])
    {
      // Populate the collision count arrays.
      for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
      {
        intraProcCollisions[collisionType] = intraBlockNumbers[collisionType].size();
        interProcCollisions[collisionType] = interBlockNumbers[collisionType].size();
      }

      // Data about local sites.
      localFluidSites = 0;

      // Data about contiguous local sites. First intra-proc stuff, then inter-proc.
      for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
      {
        for (unsigned indexInType = 0; indexInType < intraProcCollisions[collisionType]; indexInType++)
        {
          siteData.push_back(intraSiteData[collisionType][indexInType]);
          wallNormalAtSite.push_back(intraWallNormals[collisionType][indexInType]);
          for (Direction direction = 1; direction < D3Q15::NUMVECTORS; direction++)
          {
            distanceToWall.push_back(intraWallDistance[collisionType][indexInType
                * (D3Q15::NUMVECTORS - 1) + direction - 1]);
          }

          site_t blockId = intraBlockNumbers[collisionType][indexInType];
          site_t siteId = intraSiteNumbers[collisionType][indexInType];

          Blocks[blockId].localContiguousIndex[siteId] = localFluidSites;
          localFluidSites++;
        }
      }

      for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
      {
        for (unsigned indexInType = 0; indexInType < interProcCollisions[collisionType]; indexInType++)
        {
          siteData.push_back(interSiteData[collisionType][indexInType]);
          wallNormalAtSite.push_back(interWallNormals[collisionType][indexInType]);
          for (Direction direction = 1; direction < D3Q15::NUMVECTORS; direction++)
          {
            distanceToWall.push_back(interWallDistance[collisionType][indexInType
                * (D3Q15::NUMVECTORS - 1) + direction - 1]);
          }

          site_t blockId = interBlockNumbers[collisionType][indexInType];
          site_t siteId = interSiteNumbers[collisionType][indexInType];

          Blocks[blockId].localContiguousIndex[siteId] = localFluidSites;
          localFluidSites++;
        }
      }

      fOld.resize(localFluidSites * D3Q15::NUMVECTORS + 1 + totalSharedFs);
      fNew.resize(localFluidSites * D3Q15::NUMVECTORS + 1 + totalSharedFs);
    }

    void LatticeData::CollectFluidSiteDistribution()
    {
      fluidSitesOnEachProcessor.resize(topology::NetworkTopology::Instance()->GetProcessorCount());

      hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Gathering lattice info.");
      MPI_Allgather(&localFluidSites,
                    1,
                    MpiDataType<site_t> (),
                    &fluidSitesOnEachProcessor[0],
                    1,
                    MpiDataType<site_t> (),
                    MPI_COMM_WORLD);

      totalFluidSites = 0;
      for (proc_t ii = 0; ii < topology::NetworkTopology::Instance()->GetProcessorCount(); ++ii)
      {
        totalFluidSites += fluidSitesOnEachProcessor[ii];
      }
    }

    void LatticeData::CollectGlobalSiteExtrema()
    {
      site_t localMins[3];
      site_t localMaxes[3];
      localMins[0] = std::numeric_limits<site_t>::max();
      localMins[1] = std::numeric_limits<site_t>::max();
      localMins[2] = std::numeric_limits<site_t>::max();
      localMaxes[0] = 0;
      localMaxes[1] = 0;
      localMaxes[2] = 0;
      for (site_t siteI = 0; siteI < GetXSiteCount(); ++siteI)
      {
        for (site_t siteJ = 0; siteJ < GetYSiteCount(); ++siteJ)
        {
          for (site_t siteK = 0; siteK < GetZSiteCount(); ++siteK)
          {
            const proc_t *procId = GetProcIdFromGlobalCoords(util::Vector3D<site_t>(siteI,
                                                                                    siteJ,
                                                                                    siteK));
            if (procId == NULL || *procId != topology::NetworkTopology::Instance()->GetLocalRank())
            {
              continue;
            }
            else
            {
              localMins[0] = hemelb::util::NumericalFunctions::min(localMins[0], siteI);
              localMins[1] = hemelb::util::NumericalFunctions::min(localMins[1], siteJ);
              localMins[2] = hemelb::util::NumericalFunctions::min(localMins[2], siteK);
              localMaxes[0] = hemelb::util::NumericalFunctions::max(localMaxes[0], siteI);
              localMaxes[1] = hemelb::util::NumericalFunctions::max(localMaxes[1], siteJ);
              localMaxes[2] = hemelb::util::NumericalFunctions::max(localMaxes[2], siteK);
            }
          }

        }

      }

      site_t siteMins[3], siteMaxes[3];
      MPI_Allreduce(localMins, siteMins, 3, MpiDataType<site_t> (), MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, siteMaxes, 3, MpiDataType<site_t> (), MPI_MAX, MPI_COMM_WORLD);
      for (unsigned ii = 0; ii < 3; ++ii)
      {
        globalSiteMins[ii] = siteMins[ii];
        globalSiteMaxes[ii] = siteMaxes[ii];
      }
    }

    void LatticeData::CleanEmptyBlocks()
    {
      std::vector<bool> blockIsOnThisRank(blockCount, false);

      for (site_t blockNumber = 0; blockNumber < GetBlockCount(); blockNumber++)
      {
        const BlockData& currentDataBlock = Blocks[blockNumber];

        // If we are in a block of solids, move to the next block.
        if (currentDataBlock.localContiguousIndex.size() == 0)
        {
          continue;
        }

        for (site_t localSiteIndex = 0; localSiteIndex < GetSitesPerBlockVolumeUnit(); localSiteIndex++)
        {
          if (topology::NetworkTopology::Instance()->GetLocalRank()
              == currentDataBlock.processorRankForEachBlockSite[localSiteIndex])
          {
            blockIsOnThisRank[blockNumber] = true;
            break;
          }
        }
      }

      // If we are in a block of solids, we set map_block[n].site_data to NULL.
      for (site_t n = 0; n < GetBlockCount(); n++)
      {
        if (!blockIsOnThisRank[n])
        {
          Blocks[n] = BlockData();
          continue;
        }
      }
    }

    void LatticeData::InitialiseNeighbourStuff()
    {
      // Allocate the index in which to put the distribution functions received from the other
      // process.
      std::vector < std::vector<site_t> > sharedFLocationForEachProc = std::vector<std::vector<
          site_t> >(topology::NetworkTopology::Instance()->GetProcessorCount());

      site_t totalSharedFsSoFar = 0;

      // Set the remaining neighbouring processor data.
      for (size_t n = 0; n < neighbouringProcs.size(); n++)
      {
        // Pointing to a few things, but not setting any variables.
        // FirstSharedF points to start of shared_fs.
        neighbouringProcs[n].FirstSharedF = GetLocalFluidSiteCount() * D3Q15::NUMVECTORS + 1
            + totalSharedFsSoFar;
        totalSharedFsSoFar += neighbouringProcs[n].SharedFCount;
      }

      InitialiseNeighbourLookup( sharedFLocationForEachProc);

      InitialisePointToPointComms(sharedFLocationForEachProc);

      InitialiseReceiveLookup(sharedFLocationForEachProc);
    }

    void LatticeData::InitialiseNeighbourLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc)
    {
      const proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();

      neighbourIndices.resize(D3Q15::NUMVECTORS * localFluidSites);

      site_t n = -1;

      // Iterate over blocks in global co-ords.
      for (site_t i = 0; i < GetXSiteCount(); i += GetBlockSize())
      {
        for (site_t j = 0; j < GetYSiteCount(); j += GetBlockSize())
        {
          for (site_t k = 0; k < GetZSiteCount(); k += GetBlockSize())
          {
            n++;
            geometry::BlockData& map_block_p = Blocks[n];
            if (map_block_p.localContiguousIndex.size() == 0)
            {
              continue;
            }
            site_t m = -1;
            // Iterate over sites within the block.
            for (site_t site_i = i; site_i < i + GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + GetBlockSize(); site_k++)
                {
                  // If a site is not on this process, continue.
                  m++;
                  if (localRank != map_block_p.processorRankForEachBlockSite[m])
                  {
                    continue;
                  }
                  // Get site data, which is the number of the fluid site on this proc..
                  site_t localIndex = map_block_p.localContiguousIndex[m];
                  // Set neighbour location for the distribution component at the centre of
                  // this site.
                  SetNeighbourLocation(localIndex, 0, localIndex * D3Q15::NUMVECTORS + 0);
                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Work out positions of neighbours.
                    site_t neighbourI = site_i + D3Q15::CX[l];
                    site_t neighbourJ = site_j + D3Q15::CY[l];
                    site_t neighbourK = site_k + D3Q15::CZ[l];
                    if (!IsValidLatticeSite(neighbourI, neighbourJ, neighbourK))
                    {
                      // Set the neighbour location to the rubbish site.
                      SetNeighbourLocation(localIndex,
                                           l,
                                           GetLocalFluidSiteCount() * D3Q15::NUMVECTORS);
                      continue;
                    }
                    // Get the id of the processor which the neighbouring site lies on.
                    const proc_t *proc_id_p =
                        GetProcIdFromGlobalCoords(util::Vector3D<site_t>(neighbourI,
                                                                         neighbourJ,
                                                                         neighbourK));
                    if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
                    {
                      // initialize f_id to the rubbish site.
                      SetNeighbourLocation(localIndex,
                                           l,
                                           GetLocalFluidSiteCount() * D3Q15::NUMVECTORS);
                      continue;
                    }
                    else
                    // If on the same proc, set f_id of the
                    // current site and direction to the
                    // site and direction that it sends to.
                    // If we check convergence, the data for
                    // each site is split into that for the
                    // current and previous cycles.
                    if (localRank == *proc_id_p)
                    {
                      // Pointer to the neighbour.
                      site_t contigSiteId = GetContiguousSiteId(neighbourI, neighbourJ, neighbourK);
                      SetNeighbourLocation(localIndex, l, contigSiteId * D3Q15::NUMVECTORS + l);
                      continue;
                    }
                    else
                    {
                      proc_t neigh_proc_index = (proc_t) (*proc_id_p);
                      // This stores some coordinates.  We
                      // still need to know the site number.
                      // neigh_proc[ n ].f_data is now
                      // set as well, since this points to
                      // f_data.  Every process has data for
                      // its neighbours which say which sites
                      // on this process are shared with the
                      // neighbour.
                      sharedFLocationForEachProc[neigh_proc_index].push_back(site_i);
                      sharedFLocationForEachProc[neigh_proc_index].push_back(site_j);
                      sharedFLocationForEachProc[neigh_proc_index].push_back(site_k);
                      sharedFLocationForEachProc[neigh_proc_index].push_back(l);
                    }

                  }
                }
              }
            }
          }
        }
      }
    }

    void LatticeData::InitialisePointToPointComms(std::vector<std::vector<site_t> >& sharedFLocationForEachProc)
    {
      proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();

      // point-to-point communications are performed to match data to be
      // sent to/receive from different partitions; in this way, the
      // communication of the locations of the interface-dependent fluid
      // sites and the identifiers of the distribution functions which
      // propagate to different partitions is avoided (only their values
      // will be communicated). It's here!
      // Allocate the request variable.
      net::Net tempNet;
      for (size_t m = 0; m < neighbouringProcs.size(); m++)
      {
        NeighbouringProcessor *neigh_proc_p = &neighbouringProcs[m];
        // One way send receive.  The lower numbered netTop->ProcessorCount send and the higher numbered ones receive.
        // It seems that, for each pair of processors, the lower numbered one ends up with its own
        // edge sites and directions stored and the higher numbered one ends up with those on the
        // other processor.
        if (neigh_proc_p->Rank > localRank)
        {
          tempNet.RequestSend(&sharedFLocationForEachProc[neigh_proc_p->Rank][0],
                              neigh_proc_p->SharedFCount * 4,
                              neigh_proc_p->Rank);
        }
        else
        {
          sharedFLocationForEachProc[neigh_proc_p->Rank].resize(neigh_proc_p->SharedFCount * 4);
          tempNet.RequestReceive(&sharedFLocationForEachProc[neigh_proc_p->Rank][0],
                                 neigh_proc_p->SharedFCount * 4,
                                 neigh_proc_p->Rank);
        }
      }

      tempNet.Send();
      tempNet.Receive();
      tempNet.Wait();
    }

    void LatticeData::InitialiseReceiveLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc)
    {
      proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();

      f_recv_iv.resize(totalSharedFs);

      site_t f_count = GetLocalFluidSiteCount() * D3Q15::NUMVECTORS;
      site_t sharedSitesSeen = 0;
      for (size_t m = 0; m < neighbouringProcs.size(); m++)
      {
        NeighbouringProcessor *neigh_proc_p = &neighbouringProcs[m];
        for (site_t n = 0; n < neigh_proc_p->SharedFCount; n++)
        {
          // Get coordinates and direction of the distribution function to be sent to another process.
          site_t *f_data_p = &sharedFLocationForEachProc[neigh_proc_p->Rank][n * 4];
          site_t i = f_data_p[0];
          site_t j = f_data_p[1];
          site_t k = f_data_p[2];
          site_t l = f_data_p[3];
          // Correct so that each process has the correct coordinates.
          if (neigh_proc_p->Rank < localRank)
          {
            i += D3Q15::CX[l];
            j += D3Q15::CY[l];
            k += D3Q15::CZ[l];
            l = D3Q15::INVERSEDIRECTIONS[l];
          }
          // Get the fluid site number of site that will send data to another process.
          site_t contigSiteId = GetContiguousSiteId(i, j, k);
          // Set f_id to the element in the send buffer that we put the updated
          // distribution functions in.
          SetNeighbourLocation(contigSiteId, (unsigned int) (l), ++f_count);
          // Set the place where we put the received distribution functions, which is
          // f_new[number of fluid site that sends, inverse direction].
          f_recv_iv[sharedSitesSeen] = contigSiteId * D3Q15::NUMVECTORS
              + D3Q15::INVERSEDIRECTIONS[l];
          ++sharedSitesSeen;
        }
      }
    }

    const util::Vector3D<distribn_t>& LatticeData::GetNormalToWall(site_t iSiteIndex) const
    {
      return wallNormalAtSite[iSiteIndex];
    }

    site_t LatticeData::GetXSiteCount() const
    {
      return sites.x;
    }

    site_t LatticeData::GetYSiteCount() const
    {
      return sites.y;
    }

    site_t LatticeData::GetZSiteCount() const
    {
      return sites.z;
    }

    site_t LatticeData::GetXBlockCount() const
    {
      return blockCounts.x;
    }

    site_t LatticeData::GetYBlockCount() const
    {
      return blockCounts.y;
    }

    site_t LatticeData::GetZBlockCount() const
    {
      return blockCounts.z;
    }

    distribn_t LatticeData::GetVoxelSize() const
    {
      return voxelSize;
    }

    const util::Vector3D<distribn_t> LatticeData::GetOrigin() const
    {
      return origin;
    }

    unsigned int LatticeData::GetLog2BlockSize() const
    {
      return log2BlockSize;
    }

    site_t LatticeData::GetBlockSize() const
    {
      return blockSize;
    }

    site_t LatticeData::GetBlockCount() const
    {
      return blockCount;
    }

    site_t LatticeData::GetSitesPerBlockVolumeUnit() const
    {
      return sitesPerBlockVolumeUnit;
    }

    site_t LatticeData::GetBlockIdFromBlockCoords(site_t i, site_t j, site_t k) const
    {
      return (i * blockCounts.y + j) * blockCounts.z + k;
    }

    const proc_t *LatticeData::GetProcIdFromGlobalCoords(const util::Vector3D<site_t> & globalSiteCoords) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      site_t blockI = globalSiteCoords.x >> log2BlockSize;
      site_t blockJ = globalSiteCoords.y >> log2BlockSize;
      site_t blockK = globalSiteCoords.z >> log2BlockSize;

      // Get the block from the block identifiers.
      const BlockData* block = &Blocks[GetBlockIdFromBlockCoords(blockI, blockJ, blockK)];

      // If an empty (solid) block is addressed, return a NULL pointer.
      if (block->processorRankForEachBlockSite.size() == 0)
      {
        return NULL;
      }
      else
      {
        // Find site coordinates within the block
        site_t siteI = globalSiteCoords.x - (blockI << log2BlockSize);
        site_t siteJ = globalSiteCoords.y - (blockJ << log2BlockSize);
        site_t siteK = globalSiteCoords.z - (blockK << log2BlockSize);

        // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
        // mProcessorsForEachBlock)
        return &block->processorRankForEachBlockSite[ ( ( (siteI << log2BlockSize) + siteJ)
            << log2BlockSize) + siteK];
      }
    }

    bool LatticeData::IsValidBlock(site_t i, site_t j, site_t k) const
    {
      return i < blockCounts.x && j < blockCounts.y && k < blockCounts.z && i > -1 && j > -1 && k
          > -1;
    }

    bool LatticeData::IsValidLatticeSite(site_t i, site_t j, site_t k) const
    {
      return i >= 0 && j >= 0 && k >= 0 && i < sites.x && j < sites.y && k < sites.z;
    }

    const BlockData* LatticeData::GetBlock(site_t blockNumber) const
    {
      return &Blocks[blockNumber];
    }

    distribn_t* LatticeData::GetFOld(site_t siteNumber)
    {
      return &fOld[siteNumber];
    }

    distribn_t *LatticeData::GetFNew(site_t siteNumber)
    {
      return &fNew[siteNumber];
    }

    const distribn_t* LatticeData::GetFOld(site_t siteNumber) const
    {
      return &fOld[siteNumber];
    }

    const distribn_t *LatticeData::GetFNew(site_t siteNumber) const
    {
      return &fNew[siteNumber];
    }

    site_t LatticeData::GetLocalFluidSiteCount() const
    {
      return localFluidSites;
    }

    SiteType LatticeData::GetSiteType(site_t iSiteIndex) const
    {
      return siteData[iSiteIndex].GetSiteType();
    }

    int LatticeData::GetBoundaryId(site_t iSiteIndex) const
    {
      return siteData[iSiteIndex].GetBoundaryId();
    }

    site_t LatticeData::GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const
    {
      return neighbourIndices[iSiteIndex * D3Q15::NUMVECTORS + iDirectionIndex];
    }

    bool LatticeData::HasBoundary(const site_t iSiteIndex, const int iDirection) const
    {
      return siteData[iSiteIndex].HasBoundary(iDirection);
    }

    double LatticeData::GetCutDistance(site_t iSiteIndex, int iDirection) const
    {
      return distanceToWall[iSiteIndex * (D3Q15::NUMVECTORS - 1) + iDirection - 1];
    }

    SiteData LatticeData::GetSiteData(site_t iSiteIndex) const
    {
      return siteData[iSiteIndex];
    }

    site_t LatticeData::GetContiguousSiteId(site_t siteI, site_t siteJ, site_t siteK) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      site_t blockI = siteI >> log2BlockSize;
      site_t blockJ = siteJ >> log2BlockSize;
      site_t blockK = siteK >> log2BlockSize;

      // Pointer to the block
      const BlockData* lBlock = &Blocks[GetBlockIdFromBlockCoords(blockI, blockJ, blockK)];

      // Find site coordinates within the block
      site_t localSiteI = siteI - (blockI << log2BlockSize);
      site_t localSiteJ = siteJ - (blockJ << log2BlockSize);
      site_t localSiteK = siteK - (blockK << log2BlockSize);

      // Return pointer to site_data[site]
      return lBlock->localContiguousIndex[ ( ( (localSiteI << log2BlockSize) + localSiteJ)
          << log2BlockSize) + localSiteK];
    }

    const util::Vector3D<site_t> LatticeData::GetGlobalCoords(site_t blockNumber,
                                                              const util::Vector3D<site_t> & localSiteCoords) const
    {
      site_t blockI, blockJ, blockK;
      GetBlockIJK(blockNumber, &blockI, &blockJ, &blockK);

      return util::Vector3D<site_t>( (blockI << log2BlockSize) + localSiteCoords.x,
                                     (blockJ << log2BlockSize) + localSiteCoords.y,
                                     (blockK << log2BlockSize) + localSiteCoords.z);
    }

    site_t LatticeData::GetInnerSiteCount() const
    {
      site_t innerSiteCount = 0;
      for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
      {
        innerSiteCount += intraProcCollisions[collisionType];
      }

      return innerSiteCount;
    }

    site_t LatticeData::GetInnerCollisionCount(unsigned int collisionType) const
    {
      return intraProcCollisions[collisionType];
    }

    site_t LatticeData::GetInterCollisionCount(unsigned int collisionType) const
    {
      return interProcCollisions[collisionType];
    }

    unsigned int LatticeData::GetCollisionType(unsigned int site_data) const
    {
      return siteData[site_data].GetCollisionType();
    }

    void LatticeData::SetNeighbourLocation(site_t iSiteIndex,
                                           unsigned int iDirection,
                                           site_t iValue)
    {
      neighbourIndices[iSiteIndex * D3Q15::NUMVECTORS + iDirection] = iValue;
    }

    void LatticeData::GetBlockIJK(site_t block, site_t* blockI, site_t* blockJ, site_t* blockK) const
    {
      *blockK = block % GetZBlockCount();
      site_t blockIJData = block / GetZBlockCount();
      *blockJ = blockIJData % GetYBlockCount();
      *blockI = blockIJData / GetYBlockCount();
    }

    void LatticeData::SendAndReceive(hemelb::net::Net *net)
    {
      for (std::vector<NeighbouringProcessor>::const_iterator it = neighbouringProcs.begin(); it
          != neighbouringProcs.end(); it++)
      {
        // Request the receive into the appropriate bit of FOld.
        net->RequestReceive<distribn_t> (GetFOld( (*it).FirstSharedF),
                                         (int) ( (*it).SharedFCount),
                                          (*it).Rank);
        // Request the send from the right bit of FNew.
        net->RequestSend<distribn_t> (GetFNew( (*it).FirstSharedF),
                                      (int) ( (*it).SharedFCount),
                                       (*it).Rank);

      }
    }

    void LatticeData::SwapOldAndNew()
    {
      fOld.swap(fNew);
    }

    void LatticeData::CopyReceived()
    {
      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      for (site_t i = 0; i < totalSharedFs; i++)
      {
        *GetFNew(f_recv_iv[i]) = *GetFOld(neighbouringProcs[0].FirstSharedF + i);
      }
    }

    const std::vector<site_t>& LatticeData::GetFluidSiteCountsOnEachProc() const
    {
      return fluidSitesOnEachProcessor;
    }

    site_t LatticeData::GetFluidSiteCountOnProc(proc_t proc) const
    {
      return fluidSitesOnEachProcessor[proc];
    }

    site_t LatticeData::GetTotalFluidSites() const
    {
      return totalFluidSites;
    }

    const util::Vector3D<site_t>& LatticeData::GetGlobalSiteMins() const
    {
      return globalSiteMins;
    }

    const util::Vector3D<site_t>& LatticeData::GetGlobalSiteMaxes() const
    {
      return globalSiteMaxes;
    }
  }
}
