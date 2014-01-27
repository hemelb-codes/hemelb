// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "geometry/decomposition/OptimisedDecomposition.h"
#include "geometry/decomposition/DecompositionWeights.h"
#include "lb/lattices/D3Q27.h"
#include "log/Logger.h"
#include "net/net.h"

namespace hemelb
{
  namespace geometry
  {
    namespace decomposition
    {
      OptimisedDecomposition::OptimisedDecomposition(
          reporting::Timers& timers, net::MpiCommunicator& comms, const Geometry& geometry,
          const lb::lattices::LatticeInfo& latticeInfo, const std::vector<proc_t>& procForEachBlock,
          const std::vector<site_t>& fluidSitesOnEachBlock) :
          timers(timers), comms(comms), geometry(geometry), latticeInfo(latticeInfo),
              procForEachBlock(procForEachBlock), fluidSitesPerBlock(fluidSitesOnEachBlock)
      {
        timers[hemelb::reporting::Timers::InitialGeometryRead].Start(); //overall dbg timing

        // Calculate the site distribution and validate if appropriate.
        PopulateSiteDistribution();

        if (ShouldValidate())
        {
          ValidateVertexDistribution();
        }

        // Calculate the firstSiteIndexOnEachBlock and validate if appropriate
        PopulateFirstSiteIndexOnEachBlock();

        if (ShouldValidate())
        {
          ValidateFirstSiteIndexOnEachBlock();
        }

        // Populate the adjacency data arrays (for ParMetis) and validate if appropriate
        idx_t localVertexCount = vtxDistribn[comms.Rank() + 1] - vtxDistribn[comms.Rank()];

        PopulateAdjacencyData(localVertexCount);

        if (ShouldValidate())
        {
          ValidateAdjacencyData(localVertexCount);
        }

        log::Logger::Log<log::Trace, log::OnePerCore>("Adj length %i", localAdjacencies.size());

        timers[hemelb::reporting::Timers::InitialGeometryRead].Stop();

        // Call parmetis.
        timers[hemelb::reporting::Timers::parmetis].Start();
        log::Logger::Log<log::Debug, log::OnePerCore>("Making the call to Parmetis");

        CallParmetis(localVertexCount);

        timers[hemelb::reporting::Timers::parmetis].Stop();
        log::Logger::Log<log::Debug, log::OnePerCore>("Parmetis has finished.");

        // Convert the ParMetis results into a nice format.
        timers[hemelb::reporting::Timers::PopulateOptimisationMovesList].Start();
        log::Logger::Log<log::Debug, log::OnePerCore>("Getting moves lists for this core.");
        PopulateMovesList();
        log::Logger::Log<log::Debug, log::OnePerCore>("Done getting moves lists for this core");
        timers[hemelb::reporting::Timers::PopulateOptimisationMovesList].Stop();
      }

      void OptimisedDecomposition::CallParmetis(idx_t localVertexCount)
      {
        // From the ParMETIS documentation:
        // --------------------------------
        // Processor Pi holds ni consecutive vertices and mi corresponding edges
        //
        // xadj[ni+1] has the cumulative number of adjacencies per vertex (with a leading 0 on each processor)
        // vwgt[ni] has vertex weight coefficients and can be NULL
        // adjncy[mi] has the adjacent vertices for each edge (using a global index, starting at 0)
        // adjwgt[mi] has edge weights and can be NULL
        // vtxdist[P+1] has an identical array of the number of the vertices on each processor, cumulatively.
        //           So Pi has vertices from vtxdist[i] to vtxdist[i+1]-1
        // wgtflag* is 0 with no weights (1 on edges, 2 on vertices, 3 on edges & vertices)
        // numflag* is 0 for C-style numbering (1 for Fortran-style)
        // ncon* is the number of weights on each vertex
        // nparts* is the number of sub-domains (partition domains) desired
        // tpwgts* is the fraction of vertex weight to apply to each sub-domain
        // ubvec* is an array of the imbalance tolerance for each vertex weight
        // options* is an int array of options
        //
        // edgecut[1] will contain the number of edges cut by the partitioning
        // part[ni] will contain the partition vector of the locally-stored vertices
        // comm* is a pointer to the MPI communicator of the processes involved

        // Initialise the partition vector.
        partitionVector = std::vector<idx_t>(localVertexCount, comms.Rank());

        // Weight all vertices evenly.
        //std::vector < idx_t > vertexWeight(localVertexCount, 1);

        // Populate the vertex weight data arrays (for ParMetis) and print out the number of different fluid sites on each core
        PopulateVertexWeightData(localVertexCount);

        // Set the weights of each partition to be even, and to sum to 1.
        idx_t desiredPartitionSize = comms.Size();

        std::vector<real_t> domainWeights(desiredPartitionSize,
                                          (real_t)(1.0) / ( (real_t)(desiredPartitionSize)));
        // A bunch of values ParMetis needs.
        idx_t noConstraints = 1;
        idx_t weightFlag = 2;
        idx_t numberingFlag = 0;
        idx_t edgesCut = 0;
        idx_t nDims = 3;
        idx_t options[4] = { 0, 0, 0, 0 };
        if (ShouldValidate())
        {
          // Specify that some options are set and that we should
          // debug everything.
          // Specify that we have set some options
          options[0] = 1;
          // From parmetis.h
          // We get timing info (1)
          // more timing info (2)
          // details of the graph-coarsening process (4)
          // info during graph refinement (8)
          // NOT info on matching (16)
          // info on communication during matching (32)
          // info on remappining (64)
          options[1] = 1 | 2 | 4 | 8 | 32 | 64;
        }
        real_t tolerance = 1.001F;
        log::Logger::Log<log::Debug, log::OnePerCore>("Calling ParMetis");
        // Reserve 1 on these vectors so that the reference to their first element
        // exists (even if it's unused).
        // Reserve on the vectors to be certain they're at least 1 in capacity (so &vector[0] works)
        partitionVector.reserve(1);
        vtxDistribn.reserve(1);
        adjacenciesPerVertex.reserve(1);
        localAdjacencies.reserve(1);
        vertexWeights.reserve(1);
        MPI_Comm communicator = comms;
        ParMETIS_V3_PartKway(&vtxDistribn[0],
         &adjacenciesPerVertex[0],
         &localAdjacencies[0],
         &vertexWeights[0],
         NULL,
         &weightFlag,
         &numberingFlag,
         &noConstraints,
         &desiredPartitionSize,
         &domainWeights[0],
         &tolerance,
         options,
         &edgesCut,
         &partitionVector[0],
         &communicator);
        /*ParMETIS_V3_PartGeomKway(&vtxDistribn[0],
                                 &adjacenciesPerVertex[0],
                                 &localAdjacencies[0],
                                 &vertexWeights[0],
                                 NULL,
                                 &weightFlag,
                                 &numberingFlag,
                                 &nDims,
                                 &vertexCoordinates[0],
                                 &noConstraints,
                                 &desiredPartitionSize,
                                 &domainWeights[0],
                                 &tolerance,
                                 options,
                                 &edgesCut,
                                 &partitionVector[0],
                                 &communicator);*/

        /** Preliminary development code to create a group communicator
         std::vector<int> localRanksInNode;
         int localBaseRank = comms.Rank() - (comms.Rank() % hemelbCoresPerNode);
         int coresInNodePartition = hemelbCoresPerNode;
         log::Logger::Log<log::Info, log::OnePerCore>("Cores per node %d.", hemelbCoresPerNode);

         if(localBaseRank + hemelbCoresPerNode > comms.Size()) {
         coresInNodePartition = comms.Size() - localBaseRank;
         }

         for(int i=0; i<coresInNodePartition; i++) {
         localRanksInNode.push_back(localBaseRank + i);
         }

         net::MpiGroup GroupIntraNode = comms.Group().Include(localRanksInNode);
         net::MpiCommunicator CommsIntraNode = comms.Create(GroupIntraNode);*/

        log::Logger::Log<log::Debug, log::OnePerCore>("ParMetis returned.");
        if (comms.Rank() == comms.Size() - 1)
        {
          log::Logger::Log<log::Info, log::OnePerCore>("ParMetis cut %d edges.", edgesCut);
          if (edgesCut < 1 && comms.Size() > 2)
          {
            throw Exception()
                << "The decomposition using ParMetis returned an edge cut of 0 even though there are multiple processes. "
                << "This means/implies that ParMETIS cannot properly decompose the system, and no properly load-balanced parallel simulation can be started.";
          }
        }
      }

      void OptimisedDecomposition::PopulateVertexWeightData(idx_t localVertexCount)
      {
        // These counters will be used later on to count the number of each type of vertex site
        int FluidSiteCounter = 0, WallSiteCounter = 0, IOSiteCounter = 0, WallIOSiteCounter = 0;
        int localweight = 1;
        //vertexWeights.push_back(0);

        //We define every architecture that we will switch on:
        std::string INTELSANDYBRIDGE = "INTELSANDYBRIDGE";
        std::string AMDBULLDOZER = "AMDBULLDOZER";
        std::string NEUTRAL = "NEUTRAL";

        // For each block (counting up by lowest site id)...
        for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; blockI++)
        {
          for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; blockJ++)
          {
            for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; blockK++)
            {
              const site_t blockNumber = geometry.GetBlockIdFromBlockCoordinates(blockI,
                                                                                 blockJ,
                                                                                 blockK);

              // Only consider sites on this processor.
              if (procForEachBlock[blockNumber] != comms.Rank())
              {
                continue;
              }

              const BlockReadResult& blockReadResult = geometry.Blocks[blockNumber];

              site_t m = -1;
              const int block_size = geometry.GetBlockSize();
              const int blockXCoord = blockI * block_size;
              const int blockYCoord = blockJ * block_size;
              const int blockZCoord = blockK * block_size;

              // ... iterate over sites within the block...
              for (site_t localSiteI = 0; localSiteI < block_size; localSiteI++)
              {
                for (site_t localSiteJ = 0; localSiteJ < block_size; localSiteJ++)
                {
                  for (site_t localSiteK = 0; localSiteK < block_size; localSiteK++)
                  {
                    ++m;

                    // ... only looking at non-solid sites...
                    if (blockReadResult.Sites[m].targetProcessor == BIG_NUMBER2)
                    {
                      continue;
                    }

                    //Getting Site ID to be able to identify site type
                    site_t localSiteId = geometry.GetSiteIdFromSiteCoordinates(localSiteI,
                                                                               localSiteJ,
                                                                               localSiteK);

                    //Switch structure which identifies site type and assigns the proper weight to each vertex

                    SiteData siteData(blockReadResult.Sites[localSiteId]);

                    switch (siteData.GetCollisionType())
                    {
                      case FLUID:
                        localweight = hemelbSiteWeights[0];
                        ++FluidSiteCounter;
                        break;

                      case WALL:
                        localweight = hemelbSiteWeights[1];
                        ++WallSiteCounter;
                        break;

                      case INLET:
                        localweight = hemelbSiteWeights[2];
                        ++IOSiteCounter;
                        break;

                      case OUTLET:
                        localweight = hemelbSiteWeights[3];
                        ++IOSiteCounter;
                        break;

                      case (INLET | WALL):
                        localweight = hemelbSiteWeights[4];
                        ++WallIOSiteCounter;
                        break;

                      case (OUTLET | WALL):
                        localweight = hemelbSiteWeights[5];
                        ++WallIOSiteCounter;
                        break;
                    }

                    vertexWeights.push_back(localweight);
                    vertexCoordinates.push_back(blockXCoord + localSiteI);
                    vertexCoordinates.push_back(blockYCoord + localSiteJ);
                    vertexCoordinates.push_back(blockZCoord + localSiteK);

                  }

                }

              }

            }

          }
        }

        int TotalCoreWeight = ( (FluidSiteCounter * hemelbSiteWeights[0])
            + (WallSiteCounter * hemelbSiteWeights[1]) + (IOSiteCounter * hemelbSiteWeights[2])
            + (WallIOSiteCounter * hemelbSiteWeights[4])) / hemelbSiteWeights[0];
        int TotalSites = FluidSiteCounter + WallSiteCounter + WallIOSiteCounter;

        log::Logger::Log<log::Debug, log::OnePerCore>("There are %u Bulk Flow Sites, %u Wall Sites, %u IO Sites, %u WallIO Sites on core %u. Total: %u (Weighted %u Points)",
                                                      FluidSiteCounter,
                                                      WallSiteCounter,
                                                      IOSiteCounter,
                                                      WallIOSiteCounter,
                                                      comms.Rank(),
                                                      TotalSites,
                                                      TotalCoreWeight);
      }

      void OptimisedDecomposition::PopulateSiteDistribution()
      {
        vtxDistribn.resize(comms.Size() + 1, 0);
        // Firstly, count the sites per processor. Do this off-by-one
        // to be compatible with ParMetis.
        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
          if (procForEachBlock[block] >= 0)
          {
            vtxDistribn[1 + procForEachBlock[block]] += (idx_t)(fluidSitesPerBlock[block]);
          }
        }

        // Now make the count cumulative, again off-by-one.
        for (proc_t rank = 0; rank < comms.Size(); ++rank)
        {
          vtxDistribn[rank + 1] += vtxDistribn[rank];
        }
      }

      void OptimisedDecomposition::PopulateFirstSiteIndexOnEachBlock()
      {
        // First calculate the lowest site index on each proc - relatively easy.
        std::vector<idx_t> firstSiteOnProc(vtxDistribn);
        // Now for each block (in ascending order), the smallest site index is the smallest site
        // index on its processor, incremented by the number of sites observed from that processor
        // so far.
        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
          proc_t proc = procForEachBlock[block];
          if (proc < 0)
          {
            firstSiteIndexPerBlock.push_back(-1);
          }
          else
          {
            firstSiteIndexPerBlock.push_back(firstSiteOnProc[proc]);
            firstSiteOnProc[proc] += (idx_t)(fluidSitesPerBlock[block]);
          }
        }

      }

      void OptimisedDecomposition::PopulateAdjacencyData(idx_t localVertexCount)
      {
        adjacenciesPerVertex.push_back(0);

        // TODO: Reimplement using traversers or iterators.
        // For each block (counting up by lowest site id)...
        for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; blockI++)
        {
          for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; blockJ++)
          {
            for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; blockK++)
            {
              const site_t blockNumber = geometry.GetBlockIdFromBlockCoordinates(blockI,
                                                                                 blockJ,
                                                                                 blockK);
              // ... considering only the ones which live on this proc...
              if (procForEachBlock[blockNumber] != comms.Rank())
              {
                continue;
              }
              const BlockReadResult& blockReadResult = geometry.Blocks[blockNumber];
              site_t m = -1;
              // ... iterate over sites within the block...
              for (site_t localSiteI = 0; localSiteI < geometry.GetBlockSize(); localSiteI++)
              {
                for (site_t localSiteJ = 0; localSiteJ < geometry.GetBlockSize(); localSiteJ++)
                {
                  for (site_t localSiteK = 0; localSiteK < geometry.GetBlockSize(); localSiteK++)
                  {
                    ++m;
                    // ... only looking at non-solid sites...
                    if (blockReadResult.Sites[m].targetProcessor == BIG_NUMBER2)
                    {
                      continue;
                    }
                    // ... for each lattice direction...
                    for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++)
                    {
                      // ... which leads to a valid neighbouring site...
                      site_t neighbourI = blockI * geometry.GetBlockSize() + localSiteI
                          + latticeInfo.GetVector(l).x;
                      site_t neighbourJ = blockJ * geometry.GetBlockSize() + localSiteJ
                          + latticeInfo.GetVector(l).y;
                      site_t neighbourK = blockK * geometry.GetBlockSize() + localSiteK
                          + latticeInfo.GetVector(l).z;
                      if (neighbourI < 0 || neighbourJ < 0 || neighbourK < 0
                          || neighbourI
                              >= (geometry.GetBlockSize() * geometry.GetBlockDimensions().x)
                          || neighbourJ
                              >= (geometry.GetBlockSize() * geometry.GetBlockDimensions().y)
                          || neighbourK
                              >= (geometry.GetBlockSize() * geometry.GetBlockDimensions().z))
                      {
                        continue;
                      }
                      // ... (that is actually being simulated and not a solid)...
                      site_t neighbourBlockI = neighbourI / geometry.GetBlockSize();
                      site_t neighbourBlockJ = neighbourJ / geometry.GetBlockSize();
                      site_t neighbourBlockK = neighbourK / geometry.GetBlockSize();
                      site_t neighbourSiteI = neighbourI % geometry.GetBlockSize();
                      site_t neighbourSiteJ = neighbourJ % geometry.GetBlockSize();
                      site_t neighbourSiteK = neighbourK % geometry.GetBlockSize();
                      site_t neighbourBlockId =
                          geometry.GetBlockIdFromBlockCoordinates(neighbourBlockI,
                                                                  neighbourBlockJ,
                                                                  neighbourBlockK);
                      const BlockReadResult& neighbourBlock = geometry.Blocks[neighbourBlockId];
                      site_t neighbourSiteId =
                          geometry.GetSiteIdFromSiteCoordinates(neighbourSiteI,
                                                                neighbourSiteJ,
                                                                neighbourSiteK);
                      if (neighbourBlock.Sites.size() == 0
                          || neighbourBlock.Sites[neighbourSiteId].targetProcessor == BIG_NUMBER2)
                      {
                        continue;
                      }
                      // Calculate the site's id over the whole geometry,
                      site_t neighGlobalSiteId = firstSiteIndexPerBlock[neighbourBlockId];
                      for (site_t neighSite = 0; neighSite < geometry.GetSitesPerBlock();
                          ++neighSite)
                      {
                        if (neighSite == neighbourSiteId)
                        {
                          break;
                        }
                        else if (neighbourBlock.Sites[neighSite].targetProcessor != BIG_NUMBER2)
                        {
                          ++neighGlobalSiteId;
                        }
                      }

                      // then add this to the list of adjacencies.
                      localAdjacencies.push_back( (idx_t)(neighGlobalSiteId));
                    }

                    // The cumulative count of adjacencies for this vertex is equal to the total
                    // number of adjacencies we've entered.
                    // NOTE: The prefix operator is correct here because
                    // the array has a leading 0 not relating to any site.
                    adjacenciesPerVertex.push_back(localAdjacencies.size());
                  }

                }

              }

            }

          }

        }

      }

      std::vector<idx_t> OptimisedDecomposition::CompileMoveData(
          std::map<site_t, site_t>& blockIdLookupByLastSiteIndex)
      {
        // Right. Let's count how many sites we're going to have to move. Count the local number of
        // sites to be moved, and collect the site id and the destination processor.
        std::vector<idx_t> moveData;
        const idx_t myLowest = vtxDistribn[comms.Rank()];
        const idx_t myHighest = vtxDistribn[comms.Rank() + 1] - 1;

        // For each local fluid site...
        for (idx_t ii = 0; ii <= (myHighest - myLowest); ++ii)
        {
          // ... if it's going elsewhere...
          if (partitionVector[ii] != comms.Rank())
          {
            // ... get its id on the local processor...
            idx_t localFluidSiteId = myLowest + ii;

            // ... find out which block it's on, using our lookup map...

            // A feature of std::map::equal_range is that if there's no equal key, both iterators
            // returned will point to the entry with the next greatest key. Since we store block
            // ids by last fluid site number, this immediately gives us the block id.
            std::pair<std::map<site_t, site_t>::iterator, std::map<site_t, site_t>::iterator> rangeMatch =
                blockIdLookupByLastSiteIndex.equal_range(localFluidSiteId);

            idx_t fluidSiteBlock = rangeMatch.first->second;

            // Check the block id is correct
            if (ShouldValidate())
            {
              if (procForEachBlock[fluidSiteBlock] < 0)
              {
                log::Logger::Log<log::Critical, log::OnePerCore>("Found block %i for site %i but this block has a processor of %i assigned",
                                                                 fluidSiteBlock,
                                                                 localFluidSiteId,
                                                                 procForEachBlock[fluidSiteBlock]);
              }
              if (firstSiteIndexPerBlock[fluidSiteBlock] > localFluidSiteId)
              {
                log::Logger::Log<log::Critical, log::OnePerCore>("Found block %i for site %i but sites on this block start at number %i",
                                                                 fluidSiteBlock,
                                                                 localFluidSiteId,
                                                                 firstSiteIndexPerBlock[fluidSiteBlock]);
              }
              if (firstSiteIndexPerBlock[fluidSiteBlock] + fluidSitesPerBlock[fluidSiteBlock] - 1
                  < localFluidSiteId)
              {
                log::Logger::Log<log::Critical, log::OnePerCore>("Found block %i for site %i but there are %i sites on this block starting at %i",
                                                                 fluidSiteBlock,
                                                                 localFluidSiteId,
                                                                 fluidSitesPerBlock[fluidSiteBlock],
                                                                 firstSiteIndexPerBlock[fluidSiteBlock]);
              }
            }

            // ... and find its site id within that block. Start by working out how many fluid sites
            // we have to pass before we arrive at the fluid site we're after...
            idx_t fluidSitesToPass = localFluidSiteId - firstSiteIndexPerBlock[fluidSiteBlock];
            idx_t siteIndex = 0;

            while (true)
            {
              // ... then keep going through the sites on the block until we've passed as many fluid
              // sites as we need to.
              if (geometry.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor != BIG_NUMBER2)
              {
                fluidSitesToPass--;
              }
              if (fluidSitesToPass < 0)
              {
                break;
              }
              siteIndex++;
            }

            // The above code could go wrong, so in debug logging mode, we do some extra tests.
            if (ShouldValidate())
            {
              // If we've ended up on an impossible block, or one that doesn't live on this rank,
              // inform the user.
              if (fluidSiteBlock >= geometry.GetBlockCount()
                  || procForEachBlock[fluidSiteBlock] != comms.Rank())
              {
                log::Logger::Log<log::Critical, log::OnePerCore>("Partition element %i wrongly assigned to block %u of %i (block on processor %i)",
                                                                 ii,
                                                                 fluidSiteBlock,
                                                                 geometry.GetBlockCount(),
                                                                 procForEachBlock[fluidSiteBlock]);
              }
              // Similarly, if we've ended up with an impossible site index, or a solid site,
              // print an error message.
              if (siteIndex >= geometry.GetSitesPerBlock()
                  || geometry.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor
                      == BIG_NUMBER2)
              {
                log::Logger::Log<log::Critical, log::OnePerCore>("Partition element %i wrongly assigned to site %u of %i (block %i%s)",
                                                                 ii,
                                                                 siteIndex,
                                                                 fluidSitesPerBlock[fluidSiteBlock],
                                                                 fluidSiteBlock,
                                                                 geometry.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor
                                                                     == BIG_NUMBER2 ?
                                                                   " and site is solid" :
                                                                   "");
              }
            }

            // Add the block, site and destination rank to our move list.
            moveData.push_back(fluidSiteBlock);
            moveData.push_back(siteIndex);
            moveData.push_back(partitionVector[ii]);
          }

        }

        return moveData;
      }

      void OptimisedDecomposition::ForceSomeBlocksOnOtherCores(
          std::vector<idx_t>& moveData,
          std::map<proc_t, std::vector<site_t> >& blockIdsIRequireFromX)
      {
        timers[hemelb::reporting::Timers::moveForcingNumbers].Start();

        net::Net netForMoveSending(comms);

        // We also need to force some data upon blocks, i.e. when they're receiving data from a new
        // block they didn't previously want to know about.
        std::map<proc_t, std::vector<site_t> > blockForcedUponX;
        std::vector<proc_t> numberOfBlocksIForceUponX(comms.Size(), 0);
        for (idx_t moveNumber = 0; moveNumber < (idx_t)(moveData.size()); moveNumber += 3)
        {
          proc_t target_proc = moveData[moveNumber + 2];
          site_t blockId = moveData[moveNumber];
          if (std::count(blockForcedUponX[target_proc].begin(),
                         blockForcedUponX[target_proc].end(),
                         blockId) == 0)
          {
            blockForcedUponX[target_proc].push_back(blockId);
            ++numberOfBlocksIForceUponX[target_proc];
            log::Logger::Log<log::Trace, log::OnePerCore>("I'm ensuring proc %i takes data about block %i",
                                                          target_proc,
                                                          blockId);
          }
        }

        // Now find how many blocks are being forced upon us from every other core.
        std::vector<proc_t> blocksForcedOnMe(comms.Size(), 0);
        log::Logger::Log<log::Debug, log::OnePerCore>("Moving forcing block numbers");
        MPI_Alltoall(&numberOfBlocksIForceUponX[0],
                     1,
                     net::MpiDataType<proc_t>(),
                     &blocksForcedOnMe[0],
                     1,
                     net::MpiDataType<proc_t>(),
                     comms);
        timers[hemelb::reporting::Timers::moveForcingNumbers].Stop();
        timers[hemelb::reporting::Timers::moveForcingData].Start();
        // Now get all the blocks being forced upon me.
        std::map<proc_t, std::vector<site_t> > blocksForcedOnMeByEachProc;
        for (proc_t otherProc = 0; otherProc < (proc_t) ( ( ( ( (comms.Size()))))); ++otherProc)
        {
          if (blocksForcedOnMe[otherProc] > 0)
          {
            blocksForcedOnMeByEachProc[otherProc] =
                std::vector<site_t>(blocksForcedOnMe[otherProc]);
            netForMoveSending.RequestReceiveV(blocksForcedOnMeByEachProc[otherProc], otherProc);
          }
          if (numberOfBlocksIForceUponX[otherProc] > 0)
          {
            netForMoveSending.RequestSendV(blockForcedUponX[otherProc], otherProc);
          }
          log::Logger::Log<log::Trace, log::OnePerCore>("I'm forcing %i blocks on proc %i.",
                                                        numberOfBlocksIForceUponX[otherProc],
                                                        otherProc);
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Moving forcing block ids");
        netForMoveSending.Dispatch();
        // Now go through every block forced upon me and add it to the list of ones I want.
        for (proc_t otherProc = 0; otherProc < (proc_t) ( ( ( ( (comms.Size()))))); ++otherProc)
        {
          if (blocksForcedOnMe[otherProc] > 0)
          {
            for (std::vector<site_t>::iterator it = blocksForcedOnMeByEachProc[otherProc].begin();
                it != blocksForcedOnMeByEachProc[otherProc].end(); ++it)
            {
              if (std::count(blockIdsIRequireFromX[otherProc].begin(),
                             blockIdsIRequireFromX[otherProc].end(),
                             *it) == 0)
              {
                blockIdsIRequireFromX[otherProc].push_back(*it);
                log::Logger::Log<log::Trace, log::OnePerCore>("I'm being forced to take block %i from proc %i",
                                                              *it,
                                                              otherProc);
              }
              // We also need to take all neighbours of the forced block from their processors.
              BlockLocation blockCoords = geometry.GetBlockCoordinatesFromBlockId(*it);
              // Iterate over every direction we might need (except 0 as we obviously already have
              // that block in the list).
              for (Direction direction = 1; direction < lb::lattices::D3Q27::NUMVECTORS;
                  ++direction)
              {
                // Calculate the putative neighbour's coordinates...
                BlockLocation neighbourCoords = blockCoords
                    + BlockLocation(lb::lattices::D3Q27::CX[direction],
                                    lb::lattices::D3Q27::CY[direction],
                                    lb::lattices::D3Q27::CZ[direction]);
                // If the neighbour is a real block...
                if (geometry.AreBlockCoordinatesValid(neighbourCoords))
                {
                  // Get the block id, and check whether it has any fluid sites...
                  site_t neighbourBlockId =
                      geometry.GetBlockIdFromBlockCoordinates(neighbourCoords.x,
                                                              neighbourCoords.y,
                                                              neighbourCoords.z);
                  proc_t neighbourBlockProc = procForEachBlock[neighbourBlockId];
                  if (neighbourBlockProc >= 0)
                  {
                    // Check whether this is a block we're already interested in from that neighbour.
                    if (std::count(blockIdsIRequireFromX[neighbourBlockProc].begin(),
                                   blockIdsIRequireFromX[neighbourBlockProc].end(),
                                   neighbourBlockId) == 0)
                    {
                      // Then add it to the list of blocks we're getting from that neighbour.
                      blockIdsIRequireFromX[neighbourBlockProc].push_back(neighbourBlockId);
                      log::Logger::Log<log::Trace, log::OnePerCore>("I need to also take block %i from proc %i",
                                                                    neighbourBlockId,
                                                                    neighbourBlockProc);
                    }
                  }

                }

              }

            }

          }

        }

        timers[hemelb::reporting::Timers::moveForcingData].Stop();
      }

      void OptimisedDecomposition::GetBlockRequirements(
          std::vector<site_t>& numberOfBlocksRequiredFrom,
          std::map<proc_t, std::vector<site_t> >& blockIdsIRequireFromX,
          std::vector<site_t>& numberOfBlocksXRequiresFromMe,
          std::map<proc_t, std::vector<site_t> >& blockIdsXRequiresFromMe)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Calculating block requirements");
        timers[hemelb::reporting::Timers::blockRequirements].Start();
        // Populate numberOfBlocksRequiredFrom
        for (proc_t otherProc = 0; otherProc < (proc_t) ( ( ( ( (comms.Size()))))); ++otherProc)
        {
          numberOfBlocksRequiredFrom[otherProc] = blockIdsIRequireFromX.count(otherProc) == 0 ?
            0 :
            blockIdsIRequireFromX[otherProc].size();
          log::Logger::Log<log::Trace, log::OnePerCore>("I require a total of %i blocks from proc %i",
                                                        numberOfBlocksRequiredFrom[otherProc],
                                                        otherProc);
        }
        // Now perform the exchange s.t. each core knows how many blocks are required of it from
        // each other core.
        MPI_Alltoall(&numberOfBlocksRequiredFrom[0],
                     1,
                     net::MpiDataType<site_t>(),
                     &numberOfBlocksXRequiresFromMe[0],
                     1,
                     net::MpiDataType<site_t>(),
                     comms);
        // Awesome. Now we need to get a list of all the blocks wanted from each core by each other
        // core.
        net::Net netForMoveSending(comms);
        for (proc_t otherProc = 0; otherProc < (proc_t) ( ( ( ( (comms.Size()))))); ++otherProc)
        {
          blockIdsXRequiresFromMe[otherProc] =
              std::vector<site_t>(numberOfBlocksXRequiresFromMe[otherProc]);
          log::Logger::Log<log::Trace, log::OnePerCore>("Proc %i requires %i blocks from me",
                                                        otherProc,
                                                        blockIdsXRequiresFromMe[otherProc].size());
          netForMoveSending.RequestReceiveV(blockIdsXRequiresFromMe[otherProc], otherProc);
          netForMoveSending.RequestSendV(blockIdsIRequireFromX[otherProc], otherProc);
        }
        netForMoveSending.Dispatch();
        timers[hemelb::reporting::Timers::blockRequirements].Stop();
      }

      void OptimisedDecomposition::ShareMoveCounts(
          std::map<site_t, idx_t>& movesForEachLocalBlock,
          std::map<proc_t, std::vector<site_t> >& blockIdsXRequiresFromMe,
          std::map<site_t, std::vector<proc_t> >& coresInterestedInEachBlock,
          std::vector<idx_t>& moveData, std::map<site_t, std::vector<idx_t> >& moveDataForEachBlock,
          std::map<proc_t, std::vector<site_t> >& blockIdsIRequireFromX,
          std::vector<idx_t>& movesForEachBlockWeCareAbout)
      {
        timers[hemelb::reporting::Timers::moveCountsSending].Start();
        // Initialise the moves for each local block to 0. This handles an edge case where a local
        // block has no moves.
        for (site_t blockId = 0; blockId < geometry.GetBlockCount(); ++blockId)
        {
          if (procForEachBlock[blockId] == comms.Rank())
          {
            movesForEachLocalBlock[blockId] = 0;
          }
        }

        for (proc_t otherProc = 0; otherProc < (proc_t) ( ( ( ( (comms.Size()))))); ++otherProc)
        {
          for (site_t blockNum = 0;
              blockNum < (site_t) ( ( ( ( (blockIdsXRequiresFromMe[otherProc].size())))));
              ++blockNum)
          {
            site_t blockId = blockIdsXRequiresFromMe[otherProc][blockNum];
            log::Logger::Log<log::Trace, log::OnePerCore>("Proc %i requires block %i from me",
                                                          otherProc,
                                                          blockId);
            if (coresInterestedInEachBlock.count(blockId) == 0)
            {
              coresInterestedInEachBlock[blockId] = std::vector<proc_t>();
            }
            coresInterestedInEachBlock[blockId].push_back(otherProc);
          }

        }

        for (site_t moveNumber = 0; moveNumber < (site_t) ( ( ( ( (moveData.size())))));
            moveNumber += 3)
        {
          site_t blockId = moveData[moveNumber];
          if (moveDataForEachBlock.count(blockId) == 0)
          {
            moveDataForEachBlock[blockId] = std::vector<idx_t>();
          }
          moveDataForEachBlock[blockId].push_back(blockId);
          moveDataForEachBlock[blockId].push_back(moveData[moveNumber + 1]);
          moveDataForEachBlock[blockId].push_back(moveData[moveNumber + 2]);
          movesForEachLocalBlock[blockId]++;
        }

        net::Net netForMoveSending(comms);
        for (proc_t otherProc = 0; otherProc < (proc_t) ( ( ( ( (comms.Size()))))); ++otherProc)
        {
          for (std::vector<site_t>::iterator it = blockIdsIRequireFromX[otherProc].begin();
              it != blockIdsIRequireFromX[otherProc].end(); ++it)
          {
            netForMoveSending.RequestReceiveR(movesForEachBlockWeCareAbout[*it], otherProc);
            log::Logger::Log<log::Trace, log::OnePerCore>("I want the move count for block %i from proc %i",
                                                          *it,
                                                          otherProc);
          }
          for (std::vector<site_t>::iterator it = blockIdsXRequiresFromMe[otherProc].begin();
              it != blockIdsXRequiresFromMe[otherProc].end(); ++it)
          {
            netForMoveSending.RequestSendR(movesForEachLocalBlock[*it], otherProc);
            log::Logger::Log<log::Trace, log::OnePerCore>("I'm sending move count for block %i to proc %i",
                                                          *it,
                                                          otherProc);
          }
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Sending move counts");
        netForMoveSending.Dispatch();
        timers[hemelb::reporting::Timers::moveCountsSending].Stop();
      }

      void OptimisedDecomposition::ShareMoveData(
          std::vector<idx_t> movesForEachBlockWeCareAbout,
          std::map<proc_t, std::vector<site_t> > blockIdsIRequireFromX,
          std::map<proc_t, std::vector<site_t> > blockIdsXRequiresFromMe,
          std::map<site_t, std::vector<idx_t> > moveDataForEachBlock)
      {
        timers[hemelb::reporting::Timers::moveDataSending].Start();
        idx_t totalMovesToReceive = 0;
        for (site_t blockId = 0; blockId < geometry.GetBlockCount(); ++blockId)
        {
          totalMovesToReceive += movesForEachBlockWeCareAbout[blockId];
        }
        log::Logger::Log<log::Trace, log::OnePerCore>("I'm expecting a total of %i moves",
                                                      totalMovesToReceive);
        // Gather the moves to the places they need to go to.
        // Moves list has block, site id, destination proc
        movesList.resize(totalMovesToReceive * 3);
        idx_t localMoveId = 0;

        net::Net netForMoveSending(comms);

        for (proc_t otherProc = 0; otherProc < (proc_t) ( ( ( ( (comms.Size()))))); ++otherProc)
        {
          allMoves[otherProc] = 0;
          for (std::vector<site_t>::iterator it = blockIdsIRequireFromX[otherProc].begin();
              it != blockIdsIRequireFromX[otherProc].end(); ++it)
          {
            if (movesForEachBlockWeCareAbout[*it] > 0)
            {
              netForMoveSending.RequestReceive(&movesList[localMoveId * 3],
                                               3 * movesForEachBlockWeCareAbout[*it],
                                               otherProc);
              localMoveId += movesForEachBlockWeCareAbout[*it];
              allMoves[otherProc] += movesForEachBlockWeCareAbout[*it];
              log::Logger::Log<log::Trace, log::OnePerCore>("Expect %i moves from from proc %i about block %i",
                                                            movesForEachBlockWeCareAbout[*it],
                                                            otherProc,
                                                            *it);
            }
          }

          for (std::vector<site_t>::iterator it = blockIdsXRequiresFromMe[otherProc].begin();
              it != blockIdsXRequiresFromMe[otherProc].end(); ++it)
          {
            if (moveDataForEachBlock[*it].size() > 0)
            {
              netForMoveSending.RequestSendV(moveDataForEachBlock[*it], otherProc);
              log::Logger::Log<log::Trace, log::OnePerCore>("Sending %i moves from to proc %i about block %i",
                                                            moveDataForEachBlock[*it].size() / 3,
                                                            otherProc,
                                                            *it);
            }
          }

          log::Logger::Log<log::Trace, log::OnePerCore>("%i moves from proc %i",
                                                        allMoves[otherProc],
                                                        otherProc);
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Sending move data");
        netForMoveSending.Dispatch();
        timers[hemelb::reporting::Timers::moveDataSending].Stop();
      }

      /**
       * Returns a list of the fluid sites to be moved.
       *
       * NOTE: This function's return value is a dynamically-allocated array of all the moves to be
       * performed, ordered by (origin processor [with a count described by the content of the first
       * parameter], site id on the origin processor). The contents of the array are contiguous
       * triplets of ints: (block id, site id on block, destination rank).
       *
       * @return
       */
      void OptimisedDecomposition::PopulateMovesList()
      {
        allMoves = std::vector<idx_t>(comms.Size());

        // Create a map for looking up block Ids: the map is from the contiguous site index
        // of the last fluid site on the block, to the block id.
        std::map<site_t, site_t> blockIdLookupByLastSiteIndex;
        for (site_t blockId = 0; blockId < geometry.GetBlockCount(); ++blockId)
        {
          if (procForEachBlock[blockId] >= 0 && procForEachBlock[blockId] != BIG_NUMBER2)
          {
            site_t lastFluidSiteId = firstSiteIndexPerBlock[blockId] + fluidSitesPerBlock[blockId]
                - 1;
            blockIdLookupByLastSiteIndex[lastFluidSiteId] = blockId;
          }
        }

        // Right. Let's count how many sites we're going to have to move. Count the local number of
        // sites to be moved, and collect the site id and the destination processor.
        std::vector<idx_t> moveData = CompileMoveData(blockIdLookupByLastSiteIndex);
        // Spread the move data around
        log::Logger::Log<log::Debug, log::OnePerCore>("Starting to spread move data");
        // First, for each core, gather a list of which blocks the current core wants to
        // know more data about.
        // Handily, the blocks we want to know about are exactly those for which we already
        // have some data.
        std::map<proc_t, std::vector<site_t> > blockIdsIRequireFromX;
        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
          if (!geometry.Blocks[block].Sites.empty())
          {
            proc_t residentProc = procForEachBlock[block];
            blockIdsIRequireFromX[residentProc].push_back(block);
            log::Logger::Log<log::Trace, log::OnePerCore>("I require block %i from proc %i (running total %i)",
                                                          block,
                                                          residentProc,
                                                          blockIdsIRequireFromX[residentProc].size());
          }
        }

        // We also need to force some data upon blocks, i.e. when they're receiving data from a new
        // block they didn't previously want to know about.
        ForceSomeBlocksOnOtherCores(moveData, blockIdsIRequireFromX);

        // Now we want to spread this info around so that each core knows which blocks each other
        // requires from it.
        std::vector<site_t> numberOfBlocksRequiredFrom(comms.Size(), 0);
        std::vector<site_t> numberOfBlocksXRequiresFromMe(comms.Size(), 0);
        std::map<proc_t, std::vector<site_t> > blockIdsXRequiresFromMe;
        GetBlockRequirements(numberOfBlocksRequiredFrom,
                             blockIdsIRequireFromX,
                             numberOfBlocksXRequiresFromMe,
                             blockIdsXRequiresFromMe);

        // OK, now to get on with the actual sending of the data...
        // Except first, it'll be helpful to organise things by blocks.
        std::map<site_t, std::vector<proc_t> > coresInterestedInEachBlock;
        std::map<site_t, std::vector<idx_t> > moveDataForEachBlock;
        std::map<site_t, idx_t> movesForEachLocalBlock;
        // And it'll also be super-handy to know how many moves we're going to have locally.
        std::vector<idx_t> movesForEachBlockWeCareAbout(geometry.GetBlockCount(), 0);
        ShareMoveCounts(movesForEachLocalBlock,
                        blockIdsXRequiresFromMe,
                        coresInterestedInEachBlock,
                        moveData,
                        moveDataForEachBlock,
                        blockIdsIRequireFromX,
                        movesForEachBlockWeCareAbout);

        ShareMoveData(movesForEachBlockWeCareAbout,
                      blockIdsIRequireFromX,
                      blockIdsXRequiresFromMe,
                      moveDataForEachBlock);
      }

      bool OptimisedDecomposition::ShouldValidate() const
      {
#ifdef HEMELB_VALIDATE_GEOMETRY
        return true;
#else
        return false;
#endif
      }

      void OptimisedDecomposition::ValidateVertexDistribution()
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the vertex distribution.");
        // vtxDistribn should be the same on all cores.
        std::vector<idx_t> vtxDistribnRecv = comms.AllReduce(vtxDistribn, MPI_MIN);

        for (proc_t rank = 0; rank < comms.Size() + 1; ++rank)
        {
          if (vtxDistribn[rank] != vtxDistribnRecv[rank])
          {
            log::Logger::Log<log::Critical, log::OnePerCore>("vtxDistribn[%i] was %li but at least one other core had it as %li.",
                                                             rank,
                                                             vtxDistribn[rank],
                                                             vtxDistribnRecv[rank]);
          }
        }

      }

      void OptimisedDecomposition::ValidateAdjacencyData(idx_t localVertexCount)
      {
        // If we're using debugging logs, check that the arguments are consistent across all cores.
        // To verify: vtxDistribn, adjacenciesPerVertex, adjacencies
        if (ShouldValidate())
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("Validating the graph adjacency structure");
          // Create an array of lists to store all of this node's adjacencies, arranged by the
          // proc the adjacent vertex is on.
          std::vector<std::multimap<idx_t, idx_t> > adjByNeighProc(comms.Size(),
                                                                   std::multimap<idx_t, idx_t>());
          // The adjacency data should correspond across all cores.
          for (idx_t index = 0; index < localVertexCount; ++index)
          {
            idx_t vertex = vtxDistribn[comms.Rank()] + index;
            // Iterate over each adjacency (of each vertex).
            for (idx_t adjNumber = 0;
                adjNumber < (adjacenciesPerVertex[index + 1] - adjacenciesPerVertex[index]);
                ++adjNumber)
            {
              idx_t adjacentVertex = localAdjacencies[adjacenciesPerVertex[index] + adjNumber];
              proc_t adjacentProc = -1;
              // Calculate the proc of the neighbouring vertex.
              for (proc_t proc = 0; proc < comms.Size(); ++proc)
              {
                if (vtxDistribn[proc] <= adjacentVertex && vtxDistribn[proc + 1] > adjacentVertex)
                {
                  adjacentProc = proc;
                  break;
                }
              }

              // If it doesn't appear to belong to any proc, something's wrong.
              if (adjacentProc == -1)
              {
                log::Logger::Log<log::Critical, log::OnePerCore>("The vertex %li has a neighbour %li which doesn\'t appear to live on any processor.",
                                                                 vertex,
                                                                 adjacentVertex);
                continue;
              }
              // Store the data if it does belong to a proc.
              adjByNeighProc[adjacentProc].insert(std::pair<idx_t, idx_t>(adjacentVertex, vertex));
            }

          }

          // Create variables for the neighbour data to go into.
          std::vector<idx_t> counts(comms.Size());
          std::vector<std::vector<idx_t> > data(comms.Size());
          log::Logger::Log<log::Debug, log::OnePerCore>("Validating neighbour data");
          // Now spread and compare the adjacency information. Larger ranks send data to smaller
          // ranks which receive the data and compare it.
          for (proc_t neigh = 0; neigh < (proc_t) ( ( ( ( (comms.Size()))))); ++neigh)
          {
            SendAdjacencyDataToLowerRankedProc(neigh,
                                               counts[neigh],
                                               data[neigh],
                                               adjByNeighProc[neigh]);
            if (neigh < comms.Rank())
            {
              // Sending arrays don't perform comparison.
              continue;
            }
            CompareAdjacencyData(neigh, counts[neigh], data[neigh], adjByNeighProc[neigh]);
          }

        }

      }

      void OptimisedDecomposition::SendAdjacencyDataToLowerRankedProc(
          proc_t neighbouringProc, idx_t& neighboursAdjacencyCount,
          std::vector<idx_t>& neighboursAdjacencyData,
          std::multimap<idx_t, idx_t>& expectedAdjacencyData)
      {
        if (neighbouringProc < comms.Rank())
        {
          // Send the array length.
          neighboursAdjacencyCount = 2 * expectedAdjacencyData.size();
          MPI_Send(&neighboursAdjacencyCount,
                   1,
                   net::MpiDataType(neighboursAdjacencyCount),
                   neighbouringProc,
                   42,
                   comms);
          // Create a sendable array (std::lists aren't organised in a sendable format).
          neighboursAdjacencyData.resize(neighboursAdjacencyCount);
          unsigned int adjacencyIndex = 0;
          for (std::multimap<idx_t, idx_t>::const_iterator it = expectedAdjacencyData.begin();
              it != expectedAdjacencyData.end(); ++it)
          {
            neighboursAdjacencyData[2 * adjacencyIndex] = it->first;
            neighboursAdjacencyData[2 * adjacencyIndex + 1] = it->second;
            ++adjacencyIndex;
          }
          // Send the data to the neighbouringProc.
          MPI_Send(&neighboursAdjacencyData[0],
                   (int) ( ( ( ( (neighboursAdjacencyCount))))),
                   net::MpiDataType<idx_t>(),
                   neighbouringProc,
                   43,
                   comms);
        }
        else
        // If this is a greater rank number than the neighbouringProc, receive the data.
        if (neighbouringProc > comms.Rank())
        {
          MPI_Recv(&neighboursAdjacencyCount,
                   1,
                   net::MpiDataType(neighboursAdjacencyCount),
                   neighbouringProc,
                   42,
                   comms,
                   MPI_STATUS_IGNORE);
          neighboursAdjacencyData.resize(neighboursAdjacencyCount);
          MPI_Recv(&neighboursAdjacencyData[0],
                   (int) ( ( ( ( (neighboursAdjacencyCount))))),
                   net::MpiDataType<idx_t>(),
                   neighbouringProc,
                   43,
                   comms,
                   MPI_STATUS_IGNORE);
        }
        else // Neigh == mTopologyRank, i.e. neighbouring vertices on the same proc
        // Duplicate the data.
        {
          neighboursAdjacencyCount = 2 * expectedAdjacencyData.size();
          neighboursAdjacencyData.resize(neighboursAdjacencyCount);
          int adjacencyIndex = 0;
          for (std::multimap<idx_t, idx_t>::const_iterator it = expectedAdjacencyData.begin();
              it != expectedAdjacencyData.end(); ++it)
          {
            neighboursAdjacencyData[2 * adjacencyIndex] = it->first;
            neighboursAdjacencyData[2 * adjacencyIndex + 1] = it->second;
            ++adjacencyIndex;
          }
        }

      }

      void OptimisedDecomposition::CompareAdjacencyData(
          proc_t neighbouringProc, idx_t neighboursAdjacencyCount,
          const std::vector<idx_t>& neighboursAdjacencyData,
          std::multimap<idx_t, idx_t>& expectedAdjacencyData)
      {
        // Now we compare. First go through the received data which is ordered as (adjacent
        // vertex, vertex) wrt the neighbouring proc.
        for (idx_t ii = 0; ii < neighboursAdjacencyCount; ii += 2)
        {
          bool found = false;
          // Go through each neighbour we know about on this proc, and check whether it
          // matches the current received neighbour-data.
          for (std::multimap<idx_t, idx_t>::iterator it =
              expectedAdjacencyData.find(neighboursAdjacencyData[ii + 1]);
              it != expectedAdjacencyData.end(); ++it)
          {
            idx_t recvAdj = it->first;
            idx_t recvAdj2 = it->second;
            if (neighboursAdjacencyData[ii] == recvAdj2
                && neighboursAdjacencyData[ii + 1] == recvAdj)
            {
              expectedAdjacencyData.erase(it);
              found = true;
              break;
            }
          }

          // No neighbour data on this proc matched the data received.
          if (!found)
          {
            log::Logger::Log<log::Critical, log::OnePerCore>("Neighbour proc %i had adjacency (%li,%li) that wasn't present on this processor.",
                                                             neighbouringProc,
                                                             neighboursAdjacencyData[ii],
                                                             neighboursAdjacencyData[ii + 1]);
          }
        }

        // The local store of adjacencies should now be empty, if there was complete matching.
        std::multimap<idx_t, idx_t>::iterator it = expectedAdjacencyData.begin();
        while (it != expectedAdjacencyData.end())
        {
          idx_t adj1 = it->first;
          idx_t adj2 = it->second;
          ++it;
          // We had neighbour-data on this proc that didn't match that received.
          log::Logger::Log<log::Critical, log::OnePerCore>("The local processor has adjacency (%li,%li) that isn't present on neighbouring processor %i.",
                                                           adj1,
                                                           adj2,
                                                           neighbouringProc);
        }
      }

      void OptimisedDecomposition::ValidateFirstSiteIndexOnEachBlock()
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the firstSiteIndexPerBlock values.");
        // Reduce finding the maximum across all nodes. Note that we have to use the maximum
        // because some cores will have -1 for a block (indicating that it has no neighbours on
        // that block.
        std::vector<idx_t> firstSiteIndexPerBlockRecv = comms.AllReduce(firstSiteIndexPerBlock,
                                                                        MPI_MAX);

        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
          if (firstSiteIndexPerBlock[block] >= 0
              && firstSiteIndexPerBlock[block] != firstSiteIndexPerBlockRecv[block])
          {
            log::Logger::Log<log::Critical, log::OnePerCore>("This core had the first site index on block %li as %li but at least one other core had it as %li.",
                                                             block,
                                                             firstSiteIndexPerBlock[block],
                                                             firstSiteIndexPerBlockRecv[block]);
          }
        }
      }
    }
  }
}
