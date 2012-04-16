#include <map>
#include <limits>

#include "debug/Debugger.h"
#include "log/Logger.h"
#include "topology/NetworkTopology.h"
#include "geometry/BlockTraverser.h"
#include "geometry/LatticeData.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    LatticeData::LatticeData(const lb::lattices::LatticeInfo& latticeInfo) :
      latticeInfo(latticeInfo)
    {
    }

    LatticeData::~LatticeData()
    {
    }

    LatticeData::LatticeData(const lb::lattices::LatticeInfo& latticeInfo, const GeometryReadResult& readResult) :
      latticeInfo(latticeInfo)
    {
      SetBasicDetails(readResult.blocks, readResult.blockSize, readResult.voxelSize, readResult.origin);

      ProcessReadSites(readResult);
      // if debugging then output beliefs regarding geometry and neighbour list
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();
        for (std::vector<NeighbouringProcessor>::iterator itNeighProc = neighbouringProcs.begin();
             itNeighProc != neighbouringProcs.end(); ++itNeighProc)
        {
          log::Logger::Log<log::Info, log::OnePerCore>(
                "LatticeData: Rank %i thinks that rank %i is a neighbour with %i shared edges\n",
                localRank, itNeighProc->Rank, itNeighProc->SharedFCount);
        }
      }
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
      blockCount = blockCounts.x * blockCounts.y * blockCounts.z;
    }

    void LatticeData::ProcessReadSites(const GeometryReadResult & readResult)
    {
      Blocks.resize(GetBlockCount());

      totalSharedFs = 0;

      std::vector<SiteData> domainEdgeSiteData[COLLISION_TYPES];
      std::vector<SiteData> midDomainSiteData[COLLISION_TYPES];
      std::vector<site_t> domainEdgeBlockNumber[COLLISION_TYPES];
      std::vector<site_t> midDomainBlockNumber[COLLISION_TYPES];
      std::vector<site_t> domainEdgeSiteNumber[COLLISION_TYPES];
      std::vector<site_t> midDomainSiteNumber[COLLISION_TYPES];
      std::vector<util::Vector3D<float> > domainEdgeWallNormals[COLLISION_TYPES];
      std::vector<util::Vector3D<float> > midDomainWallNormals[COLLISION_TYPES];
      std::vector<float> domainEdgeWallDistance[COLLISION_TYPES];
      std::vector<float> midDomainWallDistance[COLLISION_TYPES];

      proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();
      // Iterate over all blocks in site units
      for (BlockTraverser blockTraverser(*this); blockTraverser.CurrentLocationValid(); blockTraverser.TraverseOne())
      {
        site_t blockId = blockTraverser.GetCurrentIndex();
        const BlockReadResult & blockReadIn = readResult.Blocks[blockId];

        if (blockReadIn.Sites.size() == 0)
        {
          continue;
        }

        // Iterate over all sites within the current block.
        for (SiteTraverser siteTraverser = blockTraverser.GetSiteTraverser(); siteTraverser.CurrentLocationValid(); siteTraverser.TraverseOne())
        {
          site_t localSiteId = siteTraverser.GetCurrentIndex();

          if (Blocks[blockId].IsEmpty())
          {
            Blocks[blockId] = Block(GetSitesPerBlockVolumeUnit());
          }

          Blocks[blockId].SetProcessorRankForSite(localSiteId, blockReadIn.Sites[localSiteId].targetProcessor);

          // If the site is not on this processor, continue.
          if (localRank != blockReadIn.Sites[localSiteId].targetProcessor)
          {
            continue;
          }
          bool lIsInnerSite = true;
          // Iterate over all direction vectors.
          for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++)
          {
            // Find the neighbour site co-ords in this direction.
            util::Vector3D<site_t> neighbourGlobalCoords = blockTraverser.GetCurrentLocation() * readResult.blockSize
                + siteTraverser.GetCurrentLocation() + util::Vector3D<site_t>(latticeInfo.GetVector(l));

            if (neighbourGlobalCoords.x < 0 || neighbourGlobalCoords.y < 0 || neighbourGlobalCoords.z < 0
                || neighbourGlobalCoords.x >= readResult.blocks.x * readResult.blockSize || neighbourGlobalCoords.y
                >= readResult.blocks.y * readResult.blockSize || neighbourGlobalCoords.z >= readResult.blocks.z
                * readResult.blockSize)
            {
              continue;
            }

            // ... (that is actually being simulated and not a solid)...
            util::Vector3D<site_t> neighbourBlock = neighbourGlobalCoords / readResult.blockSize;
            util::Vector3D<site_t> neighbourSite = neighbourGlobalCoords % readResult.blockSize;
            site_t neighbourBlockId = readResult.GetBlockIdFromBlockCoordinates(neighbourBlock.x,
                                                                                neighbourBlock.y,
                                                                                neighbourBlock.z);

            // Move on if the neighbour is in a block of solids
            // in which case the block will contain zero sites
            // Or on if the neighbour site is solid
            // in which case the targetProcessor is BIG_NUMBER2
            // Or the neighbour is also on this processor
            // in which case the targetProcessor is localRank
            if (readResult.Blocks[neighbourBlockId].Sites.size() == 0)
            {
              continue;
            }

            site_t neighbourSiteId = readResult.GetSiteIdFromSiteCoordinates(neighbourSite.x,
                                                                             neighbourSite.y,
                                                                             neighbourSite.z);
            proc_t neighbourProc = readResult.Blocks[neighbourBlockId].Sites[neighbourSiteId].targetProcessor;
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
            proc_t lNeighbouringProcs = (proc_t) ( (neighbouringProcs.size()));
            for (proc_t mm = 0; mm < lNeighbouringProcs && flag; mm++)
            {
              // Check whether the rank for a particular neighbour has already been
              // used for this processor.  If it has, set flag to zero.
              NeighbouringProcessor* neigh_proc_p = &neighbouringProcs[mm];
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

              // if debugging then output decisions with reasoning for all neighbour processors
              if (log::Logger::ShouldDisplay<log::Debug>())
                log::Logger::Log<log::Info, log::OnePerCore>(
                     "LatticeData: added %i as neighbour for %i because site %i in block %i is neighbour to site %i in block %i in direction (%i,%i,%i)\n",
                     (int)neighbourProc, (int)localRank,
                     (int)neighbourSiteId, (int)neighbourBlockId,
                     (int)localSiteId, (int)blockId,
                     latticeInfo.GetVector(l).x,
                     latticeInfo.GetVector(l).y,
                     latticeInfo.GetVector(l).z);
            }
          }

          // Set the collision type data. map_block site data is renumbered according to
          // fluid site numbers within a particular collision type.
          SiteData siteData(blockReadIn.Sites[localSiteId]);
          int l = -1;
          switch (siteData.GetCollisionType())
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
          // We approximate the wall normal at a point by iterating through each wall intersection
          // link and seeing which one is closest.
          util::Vector3D<float> normal(NO_VALUE);
          // NB. initialised to absurdly large distance (normal scale is (0, root(3)))
          float shortestDistance = 100.0;
          for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); ++direction)
          {
            if (blockReadIn.Sites[localSiteId].links[direction - 1].type == LinkReadResult::WALL_INTERSECTION)
            {
              util::Vector3D<float> distance = latticeInfo.GetVector(direction)
                  * blockReadIn.Sites[localSiteId].links[direction - 1].distanceToIntersection;
              if (distance.GetMagnitude() < shortestDistance)
              {
                shortestDistance = distance.GetMagnitude();
                normal = distance.Normalise();
              }
            }

          }

          if (lIsInnerSite)
          {
            midDomainBlockNumber[l].push_back(blockId);
            midDomainSiteNumber[l].push_back(localSiteId);
            midDomainSiteData[l].push_back(siteData);
            midDomainWallNormals[l].push_back(normal);
            for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); direction++)
            {
              midDomainWallDistance[l].push_back(blockReadIn.Sites[localSiteId].links[direction - 1].distanceToIntersection);
            }
          }
          else
          {
            domainEdgeBlockNumber[l].push_back(blockId);
            domainEdgeSiteNumber[l].push_back(localSiteId);
            domainEdgeSiteData[l].push_back(siteData);
            domainEdgeWallNormals[l].push_back(normal);
            for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); direction++)
            {
              domainEdgeWallDistance[l].push_back(blockReadIn.Sites[localSiteId].links[direction - 1].distanceToIntersection);
            }
          }

        }

      }

      PopulateWithReadData(midDomainBlockNumber,
                           midDomainSiteNumber,
                           midDomainSiteData,
                           midDomainWallNormals,
                           midDomainWallDistance,
                           domainEdgeBlockNumber,
                           domainEdgeSiteNumber,
                           domainEdgeSiteData,
                           domainEdgeWallNormals,
                           domainEdgeWallDistance);
    }

    void LatticeData::CollectFluidSiteDistribution()
    {
      fluidSitesOnEachProcessor.resize(topology::NetworkTopology::Instance()->GetProcessorCount());
      hemelb::log::Logger::Log < hemelb::log::Warning, hemelb::log::Singleton > ("Gathering lattice info.");
      MPI_Allgather(&localFluidSites,
                    1,
                    MpiDataType<site_t>(),
                    &fluidSitesOnEachProcessor[0],
                    1,
                    MpiDataType<site_t>(),
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
      localMins[0] = std::numeric_limits < site_t > ::max();
      localMins[1] = std::numeric_limits < site_t > ::max();
      localMins[2] = std::numeric_limits < site_t > ::max();
      localMaxes[0] = 0;
      localMaxes[1] = 0;
      localMaxes[2] = 0;
      for (geometry::BlockTraverser blockSet(*this); blockSet.CurrentLocationValid(); blockSet.TraverseOne())
      {
        const geometry::Block& block = blockSet.GetCurrentBlockData();
        if (block.IsEmpty())
        {
          continue;
        }
        for (geometry::SiteTraverser siteSet = blockSet.GetSiteTraverser(); siteSet.CurrentLocationValid();
            siteSet.TraverseOne())
        {
          if (block.GetProcessorRankForSite(siteSet.GetCurrentIndex())
              == topology::NetworkTopology::Instance()->GetLocalRank())
          {
            util::Vector3D < site_t > globalCoords = blockSet.GetCurrentLocation() * GetBlockSize()
                + siteSet.GetCurrentLocation();
            localMins[0] = hemelb::util::NumericalFunctions::min(localMins[0], globalCoords.x);
            localMins[1] = hemelb::util::NumericalFunctions::min(localMins[1], globalCoords.y);
            localMins[2] = hemelb::util::NumericalFunctions::min(localMins[2], globalCoords.z);
            localMaxes[0] = hemelb::util::NumericalFunctions::max(localMaxes[0], globalCoords.x);
            localMaxes[1] = hemelb::util::NumericalFunctions::max(localMaxes[1], globalCoords.y);
            localMaxes[2] = hemelb::util::NumericalFunctions::max(localMaxes[2], globalCoords.z);
          }
        }

      }

      site_t siteMins[3], siteMaxes[3];
      MPI_Allreduce(localMins, siteMins, 3, MpiDataType<site_t>(), MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, siteMaxes, 3, MpiDataType<site_t>(), MPI_MAX, MPI_COMM_WORLD);
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
        const Block& currentDataBlock = GetBlock(blockNumber);
        // If we are in a block of solids, move to the next block.
        if (currentDataBlock.IsEmpty())
        {
          continue;
        }
        for (site_t localSiteIndex = 0; localSiteIndex < GetSitesPerBlockVolumeUnit(); localSiteIndex++)
        {
          if (topology::NetworkTopology::Instance()->GetLocalRank()
              == currentDataBlock.GetProcessorRankForSite(localSiteIndex))
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
          Blocks[n] = Block();
          continue;
        }
      }

    }

    void LatticeData::InitialiseNeighbourStuff()
    {
      // Allocate the index in which to put the distribution functions received from the other
      // process.
      std::vector < std::vector<site_t> > sharedFLocationForEachProc = std::vector < std::vector<site_t>
          > (topology::NetworkTopology::Instance()->GetProcessorCount());
      site_t totalSharedFsSoFar = 0;
      // Set the remaining neighbouring processor data.
      for (size_t n = 0; n < neighbouringProcs.size(); n++)
      {
        // Pointing to a few things, but not setting any variables.
        // FirstSharedF points to start of shared_fs.
        neighbouringProcs[n].FirstSharedF = GetLocalFluidSiteCount() * latticeInfo.GetNumVectors() + 1
            + totalSharedFsSoFar;
        totalSharedFsSoFar += neighbouringProcs[n].SharedFCount;
      }
      InitialiseNeighbourLookup (sharedFLocationForEachProc);
      InitialisePointToPointComms(sharedFLocationForEachProc);
      InitialiseReceiveLookup(sharedFLocationForEachProc);
    }

    void LatticeData::InitialiseNeighbourLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc)
    {
      const proc_t localRank = topology::NetworkTopology::Instance()->GetLocalRank();
      neighbourIndices.resize(latticeInfo.GetNumVectors() * localFluidSites);
      for (BlockTraverser blockTraverser(*this); blockTraverser.CurrentLocationValid(); blockTraverser.TraverseOne())
      {
        const Block& map_block_p = blockTraverser.GetCurrentBlockData();
        if (map_block_p.IsEmpty())
        {
          continue;
        }
        for (SiteTraverser siteTraverser = blockTraverser.GetSiteTraverser(); siteTraverser.CurrentLocationValid();
            siteTraverser.TraverseOne())
        {
          if (localRank != map_block_p.GetProcessorRankForSite(siteTraverser.GetCurrentIndex()))
          {
            continue;
          }
          // Get site data, which is the number of the fluid site on this proc..
          site_t localIndex = map_block_p.GetLocalContiguousIndexForSite(siteTraverser.GetCurrentIndex());
          // Set neighbour location for the distribution component at the centre of
          // this site.
          SetNeighbourLocation(localIndex, 0, localIndex * latticeInfo.GetNumVectors() + 0);
          for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++)
          {
            util::Vector3D < site_t > currentLocationCoords = blockTraverser.GetCurrentLocation() * blockSize
                + siteTraverser.GetCurrentLocation();
            // Work out positions of neighbours.
            util::Vector3D < site_t > neighbourCoords = currentLocationCoords + util::Vector3D < site_t
                > (latticeInfo.GetVector(l));
            if (!IsValidLatticeSite(neighbourCoords))
            {
              // Set the neighbour location to the rubbish site.
              SetNeighbourLocation(localIndex, l, GetLocalFluidSiteCount() * latticeInfo.GetNumVectors());
              continue;
            }
            // Get the id of the processor which the neighbouring site lies on.
            const proc_t proc_id_p = GetProcIdFromGlobalCoords(neighbourCoords);
            if (proc_id_p == BIG_NUMBER2)
            {
              // initialize f_id to the rubbish site.
              SetNeighbourLocation(localIndex, l, GetLocalFluidSiteCount() * latticeInfo.GetNumVectors());
              continue;
            }
            else
            // If on the same proc, set f_id of the
            // current site and direction to the
            // site and direction that it sends to.
            // If we check convergence, the data for
            // each site is split into that for the
            // current and previous cycles.
            if (localRank == proc_id_p)
            {
              // Pointer to the neighbour.
              site_t contigSiteId = GetContiguousSiteId(neighbourCoords);
              SetNeighbourLocation(localIndex, l, contigSiteId * latticeInfo.GetNumVectors() + l);
              continue;
            }
            else
            {
              // This stores some coordinates.  We
              // still need to know the site number.
              // neigh_proc[ n ].f_data is now
              // set as well, since this points to
              // f_data.  Every process has data for
              // its neighbours which say which sites
              // on this process are shared with the
              // neighbour.
              sharedFLocationForEachProc[proc_id_p].push_back(currentLocationCoords.x);
              sharedFLocationForEachProc[proc_id_p].push_back(currentLocationCoords.y);
              sharedFLocationForEachProc[proc_id_p].push_back(currentLocationCoords.z);
              sharedFLocationForEachProc[proc_id_p].push_back(l);
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
        NeighbouringProcessor* neigh_proc_p = &neighbouringProcs[m];
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
      site_t f_count = GetLocalFluidSiteCount() * latticeInfo.GetNumVectors();
      site_t sharedSitesSeen = 0;
      for (size_t m = 0; m < neighbouringProcs.size(); m++)
      {
        NeighbouringProcessor* neigh_proc_p = &neighbouringProcs[m];
        for (site_t n = 0; n < neigh_proc_p->SharedFCount; n++)
        {
          // Get coordinates and direction of the distribution function to be sent to another process.
          site_t* f_data_p = &sharedFLocationForEachProc[neigh_proc_p->Rank][n * 4];
          site_t i = f_data_p[0];
          site_t j = f_data_p[1];
          site_t k = f_data_p[2];
          site_t l = f_data_p[3];
          // Correct so that each process has the correct coordinates.
          if (neigh_proc_p->Rank < localRank)
          {
            i += latticeInfo.GetVector(l).x;
            j += latticeInfo.GetVector(l).y;
            k += latticeInfo.GetVector(l).z;
            l = latticeInfo.GetInverseIndex(l);
          }
          // Get the fluid site number of site that will send data to another process.
          util::Vector3D < site_t > location(i, j, k);
          site_t contigSiteId = GetContiguousSiteId(location);
          // Set f_id to the element in the send buffer that we put the updated
          // distribution functions in.
          SetNeighbourLocation(contigSiteId, (unsigned int) ( (l)), ++f_count);
          // Set the place where we put the received distribution functions, which is
          // f_new[number of fluid site that sends, inverse direction].
          f_recv_iv[sharedSitesSeen] = contigSiteId * latticeInfo.GetNumVectors() + latticeInfo.GetInverseIndex(l);
          ++sharedSitesSeen;
        }

      }

    }

    const util::Vector3D<distribn_t>& LatticeData::GetNormalToWall(site_t iSiteIndex) const
    {
      return wallNormalAtSite[iSiteIndex];
    }

    const lb::lattices::LatticeInfo& LatticeData::GetLatticeInfo() const
    {
      return latticeInfo;
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

    site_t LatticeData::GetLocalSiteIdFromLocalSiteCoords(const util::Vector3D<site_t>& siteCoords) const
    {
      return ( (siteCoords.x * blockSize) + siteCoords.y) * blockSize + siteCoords.z;
    }

    site_t LatticeData::GetBlockIdFromBlockCoords(const util::Vector3D<site_t>& blockCoords) const
    {
      return (blockCoords.x * blockCounts.y + blockCoords.y) * blockCounts.z + blockCoords.z;
    }

    proc_t LatticeData::GetProcIdFromGlobalCoords(const util::Vector3D<site_t>& globalSiteCoords) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      util::Vector3D<site_t> blockCoords, localSiteCoords;
      GetBlockAndLocalSiteCoords(globalSiteCoords, blockCoords, localSiteCoords);
      // Get the block from the block identifiers.
      const Block& block = GetBlock(GetBlockIdFromBlockCoords(blockCoords));
      // If an empty (solid) block is addressed, return a NULL pointer.
      if (block.IsEmpty())
      {
        return BIG_NUMBER2;
      }
      else
      {
        // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
        // mProcessorsForEachBlock)
        return block.GetProcessorRankForSite(GetLocalSiteIdFromLocalSiteCoords(localSiteCoords));
      }
    }

    bool LatticeData::IsValidBlock(site_t i, site_t j, site_t k) const
    {
      return i < blockCounts.x && j < blockCounts.y && k < blockCounts.z && i > -1 && j > -1 && k > -1;
    }

    bool LatticeData::IsValidLatticeSite(const util::Vector3D<site_t>& siteCoords) const
    {
      return siteCoords.x >= 0 && siteCoords.y >= 0 && siteCoords.z >= 0 && siteCoords.x < sites.x
          && siteCoords.y < sites.y && siteCoords.z < sites.z;
    }

    const Block& LatticeData::GetBlock(site_t blockNumber) const
    {
      return Blocks[blockNumber];
    }

    const distribn_t* LatticeData::GetFOld(site_t siteNumber) const
    {
      return &fOld[siteNumber];
    }

    distribn_t* LatticeData::GetFOld(site_t siteNumber)
    {
      return &fOld[siteNumber];
    }

    distribn_t* LatticeData::GetFNew(site_t siteNumber)
    {
      return &fNew[siteNumber];
    }

    const distribn_t* LatticeData::GetFNew(site_t siteNumber) const
    {
      return &fNew[siteNumber];
    }

    site_t LatticeData::GetLocalFluidSiteCount() const
    {
      return localFluidSites;
    }

    SiteData LatticeData::GetSiteData(site_t iSiteIndex) const
    {
      return siteData[iSiteIndex];
    }

    const util::Vector3D<site_t>& LatticeData::GetGlobalSiteCoords(site_t siteIndex) const
    {
      return globalSiteCoords[siteIndex];
    }

    site_t LatticeData::GetContiguousSiteId(util::Vector3D<site_t> location) const
    {
      // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
      util::Vector3D<site_t> blockCoords, localSiteCoords;
      GetBlockAndLocalSiteCoords(location, blockCoords, localSiteCoords);
      // Pointer to the block
      const Block& lBlock = GetBlock(GetBlockIdFromBlockCoords(blockCoords));
      // Return pointer to site_data[site]
      return lBlock.GetLocalContiguousIndexForSite(GetLocalSiteIdFromLocalSiteCoords(localSiteCoords));
    }

    const util::Vector3D<site_t> LatticeData::GetGlobalCoords(const util::Vector3D<site_t>& blockCoords,
                                                              const util::Vector3D<site_t>& localSiteCoords) const
    {
      return blockCoords * blockSize + localSiteCoords;
    }

    const util::Vector3D<site_t> LatticeData::GetGlobalCoords(site_t blockNumber,
                                                              const util::Vector3D<site_t>& localSiteCoords) const
    {
      util::Vector3D < site_t > blockCoords;
      GetBlockIJK(blockNumber, blockCoords);
      return GetGlobalCoords(blockCoords, localSiteCoords);
    }

    util::Vector3D<site_t> LatticeData::GetSiteCoordsFromSiteId(site_t siteId) const
    {
      util::Vector3D < site_t > siteCoords;
      siteCoords.z = siteId % blockSize;
      site_t siteIJData = siteId / blockSize;
      siteCoords.y = siteIJData % blockSize;
      siteCoords.x = siteIJData / blockSize;
      return siteCoords;
    }

    void LatticeData::GetBlockAndLocalSiteCoords(const util::Vector3D<site_t>& location,
                                                 util::Vector3D<site_t>& blockCoords,
                                                 util::Vector3D<site_t>& siteCoords) const
    {
      blockCoords = location / blockSize;
      siteCoords = location % blockSize;
    }

    site_t LatticeData::GetMidDomainSiteCount() const
    {
      site_t midDomainSiteCount = 0;
      for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
      {
        midDomainSiteCount += midDomainProcCollisions[collisionType];
      }
      return midDomainSiteCount;
    }

    site_t LatticeData::GetMidDomainCollisionCount(unsigned int collisionType) const
    {
      return midDomainProcCollisions[collisionType];
    }

    site_t LatticeData::GetDomainEdgeCollisionCount(unsigned int collisionType) const
    {
      return domainEdgeProcCollisions[collisionType];
    }

    void LatticeData::SetNeighbourLocation(site_t iSiteIndex, unsigned int iDirection, site_t iValue)
    {
      neighbourIndices[iSiteIndex * latticeInfo.GetNumVectors() + iDirection] = iValue;
    }

    void LatticeData::GetBlockIJK(site_t block, util::Vector3D<site_t>& blockCoords) const
    {
      blockCoords.z = block % GetZBlockCount();
      site_t blockIJData = block / GetZBlockCount();
      blockCoords.y = blockIJData % GetYBlockCount();
      blockCoords.x = blockIJData / GetYBlockCount();
    }

    void LatticeData::SendAndReceive(hemelb::net::Net* net)
    {
      for (std::vector<NeighbouringProcessor>::const_iterator it = neighbouringProcs.begin();
          it != neighbouringProcs.end(); ++it)
      {
        // Request the receive into the appropriate bit of FOld.
        net->RequestReceive < distribn_t > (GetFOld( (*it).FirstSharedF), (int) ( ( (*it).SharedFCount)), (*it).Rank);
        // Request the send from the right bit of FNew.
        net->RequestSend < distribn_t > (GetFNew( (*it).FirstSharedF), (int) ( ( (*it).SharedFCount)), (*it).Rank);

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

    Site LatticeData::GetSite(site_t localIndex)
    {
      return Site(localIndex, *this);
    }

    ConstSite LatticeData::GetSite(site_t localIndex) const
    {
      return ConstSite(localIndex, *this);
    }

    site_t LatticeData::GetTotalFluidSites() const
    {
      return totalFluidSites;
    }

    site_t LatticeData::GetFluidSiteCountOnProc(proc_t proc) const
    {
      return fluidSitesOnEachProcessor[proc];
    }

    const util::Vector3D<site_t>& LatticeData::GetGlobalSiteMins() const
    {
      return globalSiteMins;
    }

    const util::Vector3D<site_t>& LatticeData::GetGlobalSiteMaxes() const
    {
      return globalSiteMaxes;
    }

    void LatticeData::Report(ctemplate::TemplateDictionary& dictionary)
    {
      dictionary.SetIntValue("SITES", GetTotalFluidSites());
      dictionary.SetIntValue("BLOCKS",blockCount);
      dictionary.SetIntValue("SITESPERBLOCK",sitesPerBlockVolumeUnit);
      for (size_t n = 0; n < fluidSitesOnEachProcessor.size(); n++)
      {
        ctemplate::TemplateDictionary *proc = dictionary.AddSectionDictionary("PROCESSOR");
        proc->SetIntValue("RANK", n);
        proc->SetIntValue("SITES", fluidSitesOnEachProcessor[n]);
      }

    }
  }
}
