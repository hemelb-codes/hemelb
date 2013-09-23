// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <map>
#include <limits>

#include "debug/Debugger.h"
#include "log/Logger.h"
#include "net/NetworkTopology.h"
#include "geometry/BlockTraverser.h"
#include "geometry/LatticeData.h"
#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {
    LatticeData::LatticeData(const lb::lattices::LatticeInfo& latticeInfo) :
        latticeInfo(latticeInfo), neighbouringData(new neighbouring::NeighbouringLatticeData(latticeInfo))
    {
    }

    LatticeData::~LatticeData()
    {
      delete neighbouringData;
    }

    LatticeData::LatticeData(const lb::lattices::LatticeInfo& latticeInfo, const Geometry& readResult) :
        latticeInfo(latticeInfo), neighbouringData(new neighbouring::NeighbouringLatticeData(latticeInfo))
    {
      SetBasicDetails(readResult.GetBlockDimensions(),
                      readResult.GetBlockSize(),
                      readResult.GetVoxelSize(),
                      readResult.GetOrigin());

      ProcessReadSites(readResult);
      // if debugging then output beliefs regarding geometry and neighbour list
      if (log::Logger::ShouldDisplay<log::Trace>())
      {
        proc_t localRank = net::NetworkTopology::Instance()->GetLocalRank();
        for (std::vector<NeighbouringProcessor>::iterator itNeighProc = neighbouringProcs.begin();
            itNeighProc != neighbouringProcs.end(); ++itNeighProc)
        {
          log::Logger::Log<log::Trace, log::OnePerCore>("LatticeData: Rank %i thinks that rank %i is a neighbour with %i shared edges\n",
                                                        localRank,
                                                        itNeighProc->Rank,
                                                        itNeighProc->SharedDistributionCount);
        }
      }
      CollectFluidSiteDistribution();
      CollectGlobalSiteExtrema();

      InitialiseNeighbourLookups();
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

    void LatticeData::ProcessReadSites(const Geometry & readResult)
    {
      blocks.resize(GetBlockCount());

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

      proc_t localRank = net::NetworkTopology::Instance()->GetLocalRank();
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
        for (SiteTraverser siteTraverser = blockTraverser.GetSiteTraverser(); siteTraverser.CurrentLocationValid();
            siteTraverser.TraverseOne())
        {
          site_t localSiteId = siteTraverser.GetCurrentIndex();

          if (blocks[blockId].IsEmpty())
          {
            blocks[blockId] = Block(GetSitesPerBlockVolumeUnit());
          }

          blocks[blockId].SetProcessorRankForSite(localSiteId, blockReadIn.Sites[localSiteId].targetProcessor);

          // If the site is not on this processor, continue.
          if (localRank != blockReadIn.Sites[localSiteId].targetProcessor)
          {
            continue;
          }
          bool isMidDomainSite = true;
          // Iterate over all non-zero direction vectors.
          for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++)
          {
            // Find the neighbour site co-ords in this direction.
            util::Vector3D<site_t> neighbourGlobalCoords = blockTraverser.GetCurrentLocation()
                * readResult.GetBlockSize() + siteTraverser.GetCurrentLocation()
                + util::Vector3D<site_t>(latticeInfo.GetVector(l));

            if (neighbourGlobalCoords.x < 0 || neighbourGlobalCoords.y < 0 || neighbourGlobalCoords.z < 0
                || neighbourGlobalCoords.x >= readResult.GetBlockDimensions().x * readResult.GetBlockSize()
                || neighbourGlobalCoords.y >= readResult.GetBlockDimensions().y * readResult.GetBlockSize()
                || neighbourGlobalCoords.z >= readResult.GetBlockDimensions().z * readResult.GetBlockSize())
            {
              continue;
            }

            // ... (that is actually being simulated and not a solid)...
            util::Vector3D<site_t> neighbourBlock = neighbourGlobalCoords / readResult.GetBlockSize();
            util::Vector3D<site_t> neighbourSite = neighbourGlobalCoords % readResult.GetBlockSize();
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
            isMidDomainSite = false;
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
                ++neigh_proc_p->SharedDistributionCount;
                break;
              }
            }

            // If flag is 1, we didn't find a neighbour-proc with the neighbour-site on it
            // so we need a new neighbouring processor.
            if (flag)
            {
              // Store rank of neighbour in >neigh_proc[neigh_procs]
              NeighbouringProcessor lNewNeighbour;
              lNewNeighbour.SharedDistributionCount = 1;
              lNewNeighbour.Rank = neighbourProc;
              neighbouringProcs.push_back(lNewNeighbour);

              // if debugging then output decisions with reasoning for all neighbour processors
              log::Logger::Log<log::Trace, log::OnePerCore>("LatticeData: added %i as neighbour for %i because site %i in block %i is neighbour to site %i in block %i in direction (%i,%i,%i)\n",
                                                            (int) neighbourProc,
                                                            (int) localRank,
                                                            (int) neighbourSiteId,
                                                            (int) neighbourBlockId,
                                                            (int) localSiteId,
                                                            (int) blockId,
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
            case WALL:
              l = 1;
              break;
            case INLET:
              l = 2;
              break;
            case OUTLET:
              l = 3;
              break;
            case (INLET | WALL):
              l = 4;
              break;
            case (OUTLET | WALL):
              l = 5;
              break;
          }

          const util::Vector3D<float>& normal = blockReadIn.Sites[localSiteId].wallNormalAvailable ?
            blockReadIn.Sites[localSiteId].wallNormal :
            util::Vector3D<float>(NO_VALUE);

          if (isMidDomainSite)
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
      fluidSitesOnEachProcessor.resize(net::NetworkTopology::Instance()->GetProcessorCount());
      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("Gathering lattice info.");
      MPI_Allgather(&localFluidSites,
                    1,
                    net::MpiDataType<site_t>(),
                    &fluidSitesOnEachProcessor[0],
                    1,
                    net::MpiDataType<site_t>(),
                    MPI_COMM_WORLD);
      totalFluidSites = 0;
      for (proc_t ii = 0; ii < net::NetworkTopology::Instance()->GetProcessorCount(); ++ii)
      {
        totalFluidSites += fluidSitesOnEachProcessor[ii];
      }
    }

    void LatticeData::CollectGlobalSiteExtrema()
    {
      site_t localMins[3];
      site_t localMaxes[3];

      for (unsigned dim = 0; dim < 3; ++dim)
      {
        localMins[dim] = std::numeric_limits<site_t>::max();
        localMaxes[dim] = 0;
      }

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
              == net::NetworkTopology::Instance()->GetLocalRank())
          {
            util::Vector3D<site_t> globalCoords = blockSet.GetCurrentLocation() * GetBlockSize()
                + siteSet.GetCurrentLocation();

            for (unsigned dim = 0; dim < 3; ++dim)
            {
              localMins[dim] = hemelb::util::NumericalFunctions::min(localMins[dim], globalCoords[dim]);
              localMaxes[dim] = hemelb::util::NumericalFunctions::max(localMaxes[dim], globalCoords[dim]);
            }
          }
        }

      }

      site_t siteMins[3], siteMaxes[3];
      MPI_Allreduce(localMins, siteMins, 3, net::MpiDataType<site_t>(), MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, siteMaxes, 3, net::MpiDataType<site_t>(), MPI_MAX, MPI_COMM_WORLD);
      for (unsigned ii = 0; ii < 3; ++ii)
      {
        globalSiteMins[ii] = siteMins[ii];
        globalSiteMaxes[ii] = siteMaxes[ii];
      }
    }

    void LatticeData::InitialiseNeighbourLookups()
    {
      // Allocate the index in which to put the distribution functions received from the other
      // process.
      std::vector<std::vector<site_t> > sharedDistributionLocationForEachProc =
          std::vector<std::vector<site_t> >(net::NetworkTopology::Instance()->GetProcessorCount());
      site_t totalSharedDistributionsSoFar = 0;
      // Set the remaining neighbouring processor data.
      for (size_t neighbourId = 0; neighbourId < neighbouringProcs.size(); neighbourId++)
      {
        // Pointing to a few things, but not setting any variables.
        // FirstSharedF points to start of shared_fs.
        neighbouringProcs[neighbourId].FirstSharedDistribution = GetLocalFluidSiteCount() * latticeInfo.GetNumVectors()
            + 1 + totalSharedDistributionsSoFar;
        totalSharedDistributionsSoFar += neighbouringProcs[neighbourId].SharedDistributionCount;
      }
      InitialiseNeighbourLookup(sharedDistributionLocationForEachProc);
      InitialisePointToPointComms(sharedDistributionLocationForEachProc);
      InitialiseReceiveLookup(sharedDistributionLocationForEachProc);
    }

    void LatticeData::InitialiseNeighbourLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc)
    {
      const proc_t localRank = net::NetworkTopology::Instance()->GetLocalRank();
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
          for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); direction++)
          {
            util::Vector3D<site_t> currentLocationCoords = blockTraverser.GetCurrentLocation() * blockSize
                + siteTraverser.GetCurrentLocation();
            // Work out positions of neighbours.
            util::Vector3D<site_t> neighbourCoords = currentLocationCoords
                + util::Vector3D<site_t>(latticeInfo.GetVector(direction));
            if (!IsValidLatticeSite(neighbourCoords))
            {
              // Set the neighbour location to the rubbish site.
              SetNeighbourLocation(localIndex, direction, GetLocalFluidSiteCount() * latticeInfo.GetNumVectors());
              continue;
            }
            // Get the id of the processor which the neighbouring site lies on.
            const proc_t proc_id_p = GetProcIdFromGlobalCoords(neighbourCoords);
            if (proc_id_p == BIG_NUMBER2)
            {
              // initialize f_id to the rubbish site.
              SetNeighbourLocation(localIndex, direction, GetLocalFluidSiteCount() * latticeInfo.GetNumVectors());
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
              SetNeighbourLocation(localIndex, direction, contigSiteId * latticeInfo.GetNumVectors() + direction);
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
              sharedFLocationForEachProc[proc_id_p].push_back(direction);
            }
          }
        }

      }

    }

    void LatticeData::InitialisePointToPointComms(std::vector<std::vector<site_t> >& sharedFLocationForEachProc)
    {
      proc_t localRank = net::NetworkTopology::Instance()->GetLocalRank();
      // point-to-point communications are performed to match data to be
      // sent to/receive from different partitions; in this way, the
      // communication of the locations of the interface-dependent fluid
      // sites and the identifiers of the distribution functions which
      // propagate to different partitions is avoided (only their values
      // will be communicated). It's here!
      // Allocate the request variable.
      net::Net tempNet;
      for (size_t neighbourId = 0; neighbourId < neighbouringProcs.size(); neighbourId++)
      {
        NeighbouringProcessor* neigh_proc_p = &neighbouringProcs[neighbourId];
        // One way send receive.  The lower numbered netTop->ProcessorCount send and the higher numbered ones receive.
        // It seems that, for each pair of processors, the lower numbered one ends up with its own
        // edge sites and directions stored and the higher numbered one ends up with those on the
        // other processor.
        if (neigh_proc_p->Rank > localRank)
        {
          tempNet.RequestSendV(sharedFLocationForEachProc[neigh_proc_p->Rank], neigh_proc_p->Rank);
        }
        else
        {
          sharedFLocationForEachProc[neigh_proc_p->Rank].resize(neigh_proc_p->SharedDistributionCount * 4);
          tempNet.RequestReceiveV(sharedFLocationForEachProc[neigh_proc_p->Rank], neigh_proc_p->Rank);
        }
      }

      tempNet.Dispatch();
    }

    void LatticeData::InitialiseReceiveLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc)
    {
      proc_t localRank = net::NetworkTopology::Instance()->GetLocalRank();
      streamingIndicesForReceivedDistributions.resize(totalSharedFs);
      site_t f_count = GetLocalFluidSiteCount() * latticeInfo.GetNumVectors();
      site_t sharedSitesSeen = 0;
      for (size_t neighbourId = 0; neighbourId < neighbouringProcs.size(); neighbourId++)
      {
        NeighbouringProcessor* neigh_proc_p = &neighbouringProcs[neighbourId];
        for (site_t sharedDistributionId = 0; sharedDistributionId < neigh_proc_p->SharedDistributionCount;
            sharedDistributionId++)
        {
          // Get coordinates and direction of the distribution function to be sent to another process.
          site_t* f_data_p = &sharedFLocationForEachProc[neigh_proc_p->Rank][sharedDistributionId * 4];
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
          util::Vector3D<site_t> location(i, j, k);
          site_t contigSiteId = GetContiguousSiteId(location);
          // Set f_id to the element in the send buffer that we put the updated
          // distribution functions in.
          SetNeighbourLocation(contigSiteId, (unsigned int) ( (l)), ++f_count);
          // Set the place where we put the received distribution functions, which is
          // f_new[number of fluid site that sends, inverse direction].
          streamingIndicesForReceivedDistributions[sharedSitesSeen] = contigSiteId * latticeInfo.GetNumVectors()
              + latticeInfo.GetInverseIndex(l);
          ++sharedSitesSeen;
        }

      }

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
      if (i < 0 || i >= blockCounts.x)
        return false;
      if (j < 0 || j >= blockCounts.y)
        return false;
      if (k < 0 || k >= blockCounts.z)
        return false;

      return true;
    }

    bool LatticeData::IsValidBlock(const util::Vector3D<site_t>& blockCoords) const
    {
      if (blockCoords.x < 0 || blockCoords.x >= blockCounts.x)
        return false;
      if (blockCoords.y < 0 || blockCoords.y >= blockCounts.y)
        return false;
      if (blockCoords.z < 0 || blockCoords.z >= blockCounts.z)
        return false;

      return true;
    }

    bool LatticeData::IsValidLatticeSite(const util::Vector3D<site_t>& siteCoords) const
    {
      if (siteCoords.x < 0 || siteCoords.x >= sites.x)
        return false;
      if (siteCoords.y < 0 || siteCoords.y >= sites.y)
        return false;
      if (siteCoords.z < 0 || siteCoords.z >= sites.z)
        return false;

      return true;
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

    bool LatticeData::GetContiguousSiteId(const util::Vector3D<site_t>& globalLocation,
                                          proc_t& procId,
                                          site_t& siteId) const
    {
      // convert global coordinates to local coordinates - i.e.
      // to location of block and location of site within block
      util::Vector3D<site_t> blockCoords, localSiteCoords;
      GetBlockAndLocalSiteCoords(globalLocation, blockCoords, localSiteCoords);
      if (!IsValidBlock(blockCoords) || !IsValidLatticeSite(localSiteCoords))
        return false;

      // get information for the block using the block location
      const Block& block = GetBlock(GetBlockIdFromBlockCoords(blockCoords));
      if (block.IsEmpty())
        return false;

      // get the local site id, i.e. its index within the block
      site_t localSiteIndex = GetLocalSiteIdFromLocalSiteCoords(localSiteCoords);

      // get the rank of the processor that owns the site
      procId = block.GetProcessorRankForSite(localSiteIndex);
      if (procId != net::NetworkTopology::Instance()->GetLocalRank())
        return false;
      if (procId == BIG_NUMBER2) // means that the site is solid
        return false;

      // we only know enough information to determine solid/fluid for local sites
      // get the local contiguous index of the fluid site
      if (block.SiteIsSolid(localSiteIndex))
        return false;
      else
        siteId = block.GetLocalContiguousIndexForSite(localSiteIndex);

      return true;
    }

    const util::Vector3D<site_t> LatticeData::GetGlobalCoords(site_t blockNumber,
                                                              const util::Vector3D<site_t>& localSiteCoords) const
    {
      util::Vector3D<site_t> blockCoords;
      GetBlockIJK(blockNumber, blockCoords);
      return GetGlobalCoords(blockCoords, localSiteCoords);
    }

    util::Vector3D<site_t> LatticeData::GetSiteCoordsFromSiteId(site_t siteId) const
    {
      util::Vector3D<site_t> siteCoords;
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

    void LatticeData::GetBlockIJK(site_t block, util::Vector3D<site_t>& blockCoords) const
    {
      blockCoords.z = block % blockCounts.z;
      site_t blockIJData = block / blockCounts.z;
      blockCoords.y = blockIJData % blockCounts.y;
      blockCoords.x = blockIJData / blockCounts.y;
    }

    void LatticeData::SendAndReceive(hemelb::net::Net* net)
    {
      for (std::vector<NeighbouringProcessor>::const_iterator it = neighbouringProcs.begin();
          it != neighbouringProcs.end(); ++it)
      {
        // Request the receive into the appropriate bit of FOld.
        net->RequestReceive<distribn_t>(GetFOld( (*it).FirstSharedDistribution),
                                        (int) ( ( (*it).SharedDistributionCount)),
                                        (*it).Rank);
        // Request the send from the right bit of FNew.
        net->RequestSend<distribn_t>(GetFNew( (*it).FirstSharedDistribution),
                                     (int) ( ( (*it).SharedDistributionCount)),
                                     (*it).Rank);

      }
    }

    void LatticeData::CopyReceived()
    {
      // Copy the distribution functions received from the neighbouring
      // processors into the destination buffer "f_new".
      for (site_t i = 0; i < totalSharedFs; i++)
      {
        *GetFNew(streamingIndicesForReceivedDistributions[i]) = *GetFOld(neighbouringProcs[0].FirstSharedDistribution
            + i);
      }
    }

    void LatticeData::Report(ctemplate::TemplateDictionary& dictionary)
    {
      dictionary.SetIntValue("SITES", GetTotalFluidSites());
      dictionary.SetIntValue("BLOCKS", blockCount);
      dictionary.SetIntValue("SITESPERBLOCK", sitesPerBlockVolumeUnit);
      for (size_t n = 0; n < fluidSitesOnEachProcessor.size(); n++)
      {
        ctemplate::TemplateDictionary *proc = dictionary.AddSectionDictionary("PROCESSOR");
        proc->SetIntValue("RANK", n);
        proc->SetIntValue("SITES", fluidSitesOnEachProcessor[n]);
      }
    }
    neighbouring::NeighbouringLatticeData &LatticeData::GetNeighbouringData()
    {
      return *neighbouringData;
    }
    neighbouring::NeighbouringLatticeData const & LatticeData::GetNeighbouringData() const
    {
      return *neighbouringData;
    }
  }
}
