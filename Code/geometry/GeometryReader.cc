#include <math.h>
#include <parmetis.h>

#include "io/XdrMemReader.h"
#include "geometry/LatticeData.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace geometry
  {

    // TODO This file is generally ugly. Integrate with the functions in Net which initialise the LatDat.
    // Once the interface to this object is nice and clean, we can tidy up the code here.
    LatticeData::GeometryReader::GeometryReader(const bool reserveSteeringCore)
    {
      // Each rank needs to know its global rank
      MPI_Comm_rank(MPI_COMM_WORLD, &mGlobalRank);

      // Get the group of all procs.
      MPI_Group lWorldGroup;
      MPI_Comm_group(MPI_COMM_WORLD, &lWorldGroup);

      mParticipateInTopology = !reserveSteeringCore || mGlobalRank != 0;

      // Create our own group, without the root node.
      if (reserveSteeringCore)
      {
        int lExclusions[1] = { 0 };
        MPI_Group_excl(lWorldGroup, 1, lExclusions, &mTopologyGroup);
      }
      else
      {
        mTopologyGroup = lWorldGroup;
      }

      // Create a communicator just for the domain decomposition.
      MPI_Comm_create(MPI_COMM_WORLD, mTopologyGroup, &mTopologyComm);

      // Each rank needs to know its rank wrt the domain
      // decomposition.
      if (mParticipateInTopology)
      {
        int temp = 0;
        MPI_Comm_rank(mTopologyComm, &mTopologyRank);
        MPI_Comm_size(mTopologyComm, &temp);
        mTopologySize = (unsigned int) temp;
      }
      else
      {
        mTopologyRank = -1;
        mTopologySize = 0;
      }
    }

    LatticeData::GeometryReader::~GeometryReader()
    {
      MPI_Group_free(&mTopologyGroup);

      // Note that on rank 0, this is the same as MPI_COMM_WORLD.
      if (mParticipateInTopology)
      {
        MPI_Comm_free(&mTopologyComm);
      }
    }

    void LatticeData::GeometryReader::LoadAndDecompose(GlobalLatticeData* bGlobLatDat,
                                                       lb::LbmParameters* bLbmParams,
                                                       SimConfig* bSimConfig,
                                                       double* oReadTime,
                                                       double* oOptimiseTime)
    {
      double lStart = util::myClock();

      MPI_File lFile;
      MPI_Info lFileInfo;
      int lError;

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      MPI_Info_create(&lFileInfo);

      // Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
      std::string accessStyle = "access_style";
      std::string accessStyleValue = "sequential";
      std::string buffering = "collective_buffering";
      std::string bufferingValue = "true";

      MPI_Info_set(lFileInfo, &accessStyle[0], &accessStyleValue[0]);
      MPI_Info_set(lFileInfo, &buffering[0], &bufferingValue[0]);

      // Open the file.
      lError = MPI_File_open(MPI_COMM_WORLD,
                             &bSimConfig->DataFilePath[0],
                             MPI_MODE_RDONLY,
                             lFileInfo,
                             &lFile);

      if (lError != 0)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Unable to open file %s, exiting",
                                                     bSimConfig->DataFilePath.c_str());
        fflush(0x0);
        exit(0x0);
      }
      else
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Opened config file %s",
                                                      bSimConfig->DataFilePath.c_str());
      }
      fflush(NULL);

      // Set the view to the file.
      std::string lMode = "native";
      MPI_File_set_view(lFile, 0, MPI_CHAR, MPI_CHAR, &lMode[0], MPI_INFO_NULL);

      // Read the file preamble.
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading file preamble");
      ReadPreamble(lFile, bLbmParams, bGlobLatDat);

      // Read the file header.
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading file header");

      site_t* sitesPerBlock = new site_t[bGlobLatDat->GetBlockCount()];
      unsigned int* bytesPerBlock = new unsigned int[bGlobLatDat->GetBlockCount()];

      ReadHeader(lFile, bGlobLatDat->GetBlockCount(), sitesPerBlock, bytesPerBlock);

      // Perform an initial decomposition, of which processor should read each block.
      log::Logger::Log<log::Debug, log::OnePerCore>("Beginning initial decomposition");

      proc_t* procForEachBlock = new proc_t[bGlobLatDat->GetBlockCount()];

      if (!mParticipateInTopology)
      {
        for (site_t ii = 0; ii < bGlobLatDat->GetBlockCount(); ++ii)
        {
          procForEachBlock[ii] = -1;
        }
      }
      else
      {
        BlockDecomposition(bGlobLatDat->GetBlockCount(),
                           bGlobLatDat,
                           sitesPerBlock,
                           procForEachBlock);
      }

      // Perform the initial read-in.
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading in my blocks");

      ReadInLocalBlocks(lFile, bytesPerBlock, procForEachBlock, mTopologyRank, bGlobLatDat);

      double lMiddle = util::myClock();

      // Move the file pointer to the start of the block section, after the preamble and header.
      MPI_File_seek(lFile,
                    PreambleBytes + GetHeaderLength(bGlobLatDat->GetBlockCount()),
                    MPI_SEEK_SET);

      // TODO This logic should be moved into OptimiseDomainDecomposition
      if (mParticipateInTopology)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Beginning domain decomposition optimisation");
        OptimiseDomainDecomposition(sitesPerBlock,
                                    bytesPerBlock,
                                    procForEachBlock,
                                    lFile,
                                    bGlobLatDat);
        log::Logger::Log<log::Debug, log::OnePerCore>("Ending domain decomposition optimisation");
      }
      else
      {
        ReadInLocalBlocks(lFile, bytesPerBlock, procForEachBlock, mTopologyRank, bGlobLatDat);
      }

      // Finish up - close the file, set the timings, deallocate memory.
      MPI_File_close(&lFile);
      MPI_Info_free(&lFileInfo);

      *oReadTime = lMiddle - lStart;
      *oOptimiseTime = util::myClock() - lMiddle;

      delete[] sitesPerBlock;
      delete[] bytesPerBlock;
      delete[] procForEachBlock;
    }

    /**
     * Read in the section at the beginning of the config file.
     *
     * // PREAMBLE
     * uint stress_type = load_uint();
     * uint blocks[3]; // blocks in x,y,z
     * for (int i=0; i<3; ++i)
     *        blocks[i] = load_uint();
     * // number of sites along 1 edge of a block
     * uint block_size = load_uint();
     *
     * uint nBlocks = blocks[0] * blocks[1] * blocks[2];
     *
     * double voxel_size = load_double();
     *
     * for (int i=0; i<3; ++i)
     *        site0WorldPosition[i] = load_double();
     */
    void LatticeData::GeometryReader::ReadPreamble(MPI_File xiFile,
                                                   hemelb::lb::LbmParameters * bParams,
                                                   GlobalLatticeData* bGlobalLatticeData)
    {
      // Read in the file preamble into a buffer.
      char lPreambleBuffer[PreambleBytes];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lPreambleBuffer, PreambleBytes, MPI_CHAR, &lStatus);

      // Create an Xdr translator based on the read-in data.
      hemelb::io::XdrReader preambleReader = hemelb::io::XdrMemReader(lPreambleBuffer,
                                                                      PreambleBytes);

      // Variables we'll read.
      unsigned int stressType, blocksX, blocksY, blocksZ, blockSize;
      double voxelSize, siteZeroWorldPosition[3];

      // Read in the values.
      preambleReader.readUnsignedInt(stressType);
      preambleReader.readUnsignedInt(blocksX);
      preambleReader.readUnsignedInt(blocksY);
      preambleReader.readUnsignedInt(blocksZ);
      preambleReader.readUnsignedInt(blockSize);
      preambleReader.readDouble(voxelSize);
      for (unsigned int i = 0; i < 3; ++i)
      {
        preambleReader.readDouble(siteZeroWorldPosition[i]);
      }

      // Pass the read in variables to the objects that need them.
      bParams->StressType = (lb::StressTypes) stressType;

      bGlobalLatticeData->SetBasicDetails(blocksX, blocksY, blocksZ, blockSize);
    }

    /**
     * Read the header section, with minimal information about each block.
     *
     * Note that the output is placed into the arrays sitesInEachBlock and
     * bytesUsedByBlockInDataFile, each of which must have iBlockCount
     * elements allocated.
     *
     * // HEADER
     * uint nSites[nBlocks];
     * uint nBytes[nBlocks];
     *
     * for (int i = 0; i < nBlocks; ++i) {
     *        nSites[i] = load_uint(); // sites in block
     *        nBytes[i] = load_uint(); // length, in bytes, of the block's record in this file
     * }
     */
    void LatticeData::GeometryReader::ReadHeader(MPI_File xiFile,
                                                 site_t iBlockCount,
                                                 site_t* sitesInEachBlock,
                                                 unsigned int* bytesUsedByBlockInDataFile)
    {
      site_t headerByteCount = GetHeaderLength(iBlockCount);
      // Allocate a buffer to read into, then do the reading.
      char* lHeaderBuffer = new char[headerByteCount];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lHeaderBuffer, (int) headerByteCount, MPI_CHAR, &lStatus);

      // Create a Xdr translation object to translate from binary
      hemelb::io::XdrReader preambleReader =
          hemelb::io::XdrMemReader(lHeaderBuffer, (unsigned int) headerByteCount);

      // Read in all the data.
      for (site_t ii = 0; ii < iBlockCount; ii++)
      {
        unsigned int sites, bytes;
        preambleReader.readUnsignedInt(sites);
        preambleReader.readUnsignedInt(bytes);

        sitesInEachBlock[ii] = sites;
        bytesUsedByBlockInDataFile[ii] = bytes;
      }

      delete[] lHeaderBuffer;
    }

    /**
     * Read in the necessary blocks from the file.
     *
     * BODY:
     * Block *block = new Block[nBlocks];
     * sitesPerBlock = block_size*block_size*block_size;
     *
     * for (int i = 0; i < blocks[0] * blocks[1] * blocks[2]; ++i) {
     *   if (nSites[i] == 0)
     *     // nothing
     *     continue;
     *
     *   block[i]->sites = new Site[sitesPerBlock];
     *   for (int j = 0; j < sitesPerBlock; ++j) {
     *     block[i]->sites[j]->config = load_uint();
     *     if (block[i]->sites[j]->IsAdjacentInletOrOutlet()) {
     *       for (k=0; k<3; ++k)
     *         block[i]->sites[j]->boundaryNormal[k] = load_double();
     *       block[i]->sites[j]->boundaryDistance = load_double();
     *     }
     *
     *     if (block[i]->sites[j]->IsAdjacentToWall()) {
     *       for (k=0; k<3; ++k)
     *         block[i]->sites[j]->wallNormal[k] = load_double();
     *       block[i]->sites[j]->wallDistance = load_double();
     *     }
     *
     *     for (k=0; k<14; ++k)
     *       block[i]->sites[j]->cutDistance[k] = load_double();
     *   }
     * }
     */
    void LatticeData::GeometryReader::ReadInLocalBlocks(MPI_File iFile,
                                                        const unsigned int* bytesPerBlock,
                                                        const proc_t* unitForEachBlock,
                                                        const proc_t localRank,
                                                        const GlobalLatticeData* iGlobLatDat)
    {
      // Create a list of which blocks to read in.
      bool* readBlock = new bool[iGlobLatDat->GetBlockCount()];

      for (site_t ii = 0; ii < iGlobLatDat->GetBlockCount(); ++ii)
      {
        readBlock[ii] = false;
      }

      // Read a block in if it has fluid sites and is to live on the current processor. Also read
      // in any neighbours with fluid sites.
      for (site_t blockI = 0; blockI < iGlobLatDat->GetXBlockCount(); ++blockI)
      {
        for (site_t blockJ = 0; blockJ < iGlobLatDat->GetYBlockCount(); ++blockJ)
        {
          for (site_t blockK = 0; blockK < iGlobLatDat->GetZBlockCount(); ++blockK)
          {
            site_t lBlockId = iGlobLatDat->GetBlockIdFromBlockCoords(blockI, blockJ, blockK);

            if ( (bytesPerBlock[lBlockId] == 0) || (unitForEachBlock[lBlockId] != localRank))
            {
              continue;
            }

            // Read in all neighbouring blocks with fluid sites.
            for (site_t neighI = util::NumericalFunctions::max<site_t>(0, blockI - 1); (neighI
                <= (blockI + 1)) && (neighI < iGlobLatDat->GetXBlockCount()); ++neighI)
            {
              for (site_t neighJ = util::NumericalFunctions::max<site_t>(0, blockJ - 1); (neighJ
                  <= (blockJ + 1)) && (neighJ < iGlobLatDat->GetYBlockCount()); ++neighJ)
              {
                for (site_t neighK = util::NumericalFunctions::max<site_t>(0, blockK - 1); (neighK
                    <= (blockK + 1)) && (neighK < iGlobLatDat->GetZBlockCount()); ++neighK)
                {
                  site_t lNeighId = iGlobLatDat->GetBlockIdFromBlockCoords(neighI, neighJ, neighK);

                  if (bytesPerBlock[lNeighId] > 0)
                  {
                    readBlock[lNeighId] = true;
                  }
                }
              }
            }
          }
        }
      }

      // Each site can have at most
      // * an unsigned int (config)
      // * 4 doubles
      // * 14 further doubles (in D3Q15)
      site_t maxBytesPerBlock = (iGlobLatDat->GetSitesPerBlockVolumeUnit()) * (4 * 1 + 8 * (4
          + D3Q15::NUMVECTORS - 1));
      const unsigned int BlocksToReadInOneGo = 10;
      char* readBuffer = new char[maxBytesPerBlock * BlocksToReadInOneGo];

      // This makes sure we do the right number of read operations.
      for (unsigned int readNum = 0; readNum
          <= (iGlobLatDat->GetBlockCount() / BlocksToReadInOneGo); ++readNum)
      {
        const site_t upperLimitBlockNumber =
            util::NumericalFunctions::min<site_t>(iGlobLatDat->GetBlockCount(), (readNum + 1)
                * BlocksToReadInOneGo);
        const site_t lowerLimitBlockNumber = readNum * BlocksToReadInOneGo;

        unsigned int bytesToRead = 0;
        for (site_t ii = lowerLimitBlockNumber; ii < upperLimitBlockNumber; ++ii)
        {
          bytesToRead += bytesPerBlock[ii];
        }

        MPI_Status lStatus;

        MPI_File_read_all(iFile, readBuffer, (int) bytesToRead, MPI_CHAR, &lStatus);

        io::XdrMemReader lReader(readBuffer, bytesToRead);

        for (BlockCounter lBlock(iGlobLatDat, lowerLimitBlockNumber); lBlock
            < upperLimitBlockNumber; lBlock++)
        {
          if (readBlock[lBlock] && bytesPerBlock[lBlock] > 0)
          {
            if (iGlobLatDat->Blocks[lBlock].site_data == NULL)
            {
              iGlobLatDat->Blocks[lBlock].site_data
                  = new unsigned int[iGlobLatDat->GetSitesPerBlockVolumeUnit()];
            }
            if (iGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite == NULL)
            {
              iGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite
                  = new proc_t[iGlobLatDat->GetSitesPerBlockVolumeUnit()];
            }

            site_t m = -1;

            for (site_t ii = 0; ii < iGlobLatDat->GetBlockSize(); ii++)
            {
              site_t site_i = lBlock.GetICoord(ii);

              for (site_t jj = 0; jj < iGlobLatDat->GetBlockSize(); jj++)
              {
                site_t site_j = lBlock.GetJCoord(jj);

                for (site_t kk = 0; kk < iGlobLatDat->GetBlockSize(); kk++)
                {
                  site_t site_k = lBlock.GetKCoord(kk);

                  ++m;

                  unsigned int *site_type = &iGlobLatDat->Blocks[lBlock].site_data[m];
                  if (!lReader.readUnsignedInt(*site_type))
                  {
                    std::cout << "Error reading site type\n";
                  }

                  if ( (*site_type & SITE_TYPE_MASK) == SOLID_TYPE)
                  {
                    iGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[m] = BIG_NUMBER2;
                    continue;
                  }
                  iGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[m] = -1;

                  if (iGlobLatDat->GetCollisionType(*site_type) != FLUID)
                  {
                    // Neither solid nor simple fluid
                    if (iGlobLatDat->Blocks[lBlock].wall_data == NULL)
                    {
                      iGlobLatDat->Blocks[lBlock].wall_data
                          = new WallData[iGlobLatDat->GetSitesPerBlockVolumeUnit()];
                    }

                    if (iGlobLatDat->GetCollisionType(*site_type) & INLET
                        || iGlobLatDat->GetCollisionType(*site_type) & OUTLET)
                    {
                      double temp;
                      // INLET or OUTLET or both.
                      // These values are the boundary normal and the boundary distance.
                      for (int l = 0; l < 3; l++)
                      {
                        if (!lReader.readDouble(temp))
                        {
                          std::cout << "Error reading boundary normals\n";
                        }
                      }

                      if (!lReader.readDouble(temp))
                      {
                        std::cout << "Error reading boundary distances\n";
                      }
                    }

                    if (iGlobLatDat->GetCollisionType(*site_type) & EDGE)
                    {
                      // EDGE bit set
                      for (int l = 0; l < 3; l++)
                      {
                        if (!lReader.readDouble(iGlobLatDat->Blocks[lBlock].wall_data[m].wall_nor[l]))
                        {
                          std::cout << "Error reading edge normal\n";
                        }
                      }

                      double temp;
                      if (!lReader.readDouble(temp))
                      {
                        std::cout << "Error reading edge distance\n";
                      }
                    }

                    for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
                    {
                      if (!lReader.readDouble(iGlobLatDat->Blocks[lBlock].wall_data[m].cut_dist[l]))
                      {
                        std::cout << "Error reading cut distances\n";
                      }
                    }
                  }
                } // kk
              } // jj
            } // ii
          }
          // Not being read.

          else
          {
            if (bytesPerBlock[lBlock] > 0)
            {
              if (!lReader.SetPosition(lReader.GetPosition() + bytesPerBlock[lBlock]))
              {
                std::cout << "Error setting stream position\n";
              }
            }
          }
        }
      }

      delete[] readBuffer;
      delete[] readBlock;
    }

    /**
     * Get an initial distribution of which block should be on which processor.
     *
     * To make this fast and remove a need to read in the site info about blocks,
     * we assume that neighbouring blocks with any fluid sites on have lattice links
     * between them.
     *
     * NOTE that the old version of this code used to cope with running on multiple machines,
     * by decomposing fluid sites over machines, then decomposing over processors within one machine.
     * To achieve this here, use parmetis's "nparts" parameters to first decompose over machines, then
     * to decompose within machines.
     *
     * @param initialProcForEachBlock Array of length iBlockCount, into which the rank number will be
     * written for each block
     * @param blockCountPerProc Array of length topology size, into which the number of blocks
     * allocated to each processor will be written.
     */
    void LatticeData::GeometryReader::BlockDecomposition(const site_t iBlockCount,
                                                         const GlobalLatticeData* iGlobLatDat,
                                                         const site_t* fluidSitePerBlock,
                                                         proc_t* initialProcForEachBlock)
    {
      site_t* blockCountPerProc = new site_t[mTopologySize];

      // Count of block sites.
      site_t lUnvisitedFluidBlockCount = 0;
      for (site_t ii = 0; ii < iBlockCount; ++ii)
      {
        if (fluidSitePerBlock[ii] != 0)
        {
          ++lUnvisitedFluidBlockCount;
        }
      }

      // Initialise site count per processor
      for (unsigned int ii = 0; ii < mTopologySize; ++ii)
      {
        blockCountPerProc[ii] = 0;
      }

      // Divide blocks between the processors.
      DivideBlocks(lUnvisitedFluidBlockCount,
                   iBlockCount,
                   mTopologySize,
                   blockCountPerProc,
                   initialProcForEachBlock,
                   fluidSitePerBlock,
                   iGlobLatDat);

      delete[] blockCountPerProc;
    }

    /**
     * Get an initial decomposition of the domain macro-blocks over processors (or some other
     * unit of computation, like an entire machine).
     *
     * NOTE: We need the global lattice data and fluid sites per block in order to try to keep
     * contiguous blocks together, and to skip blocks with no fluid sites.
     *
     * @param unassignedBlocks
     * @param blockCount
     * @param unitCount
     * @param blocksOnEachUnit
     * @param unitForEachBlock
     * @param fluidSitesPerBlock
     * @param iGlobLatDat
     */
    void LatticeData::GeometryReader::DivideBlocks(site_t unassignedBlocks,
                                                   site_t blockCount,
                                                   proc_t unitCount,
                                                   site_t* blocksOnEachUnit,
                                                   proc_t* unitForEachBlock,
                                                   const site_t* fluidSitesPerBlock,
                                                   const GlobalLatticeData* iGlobLatDat)
    {
      // Initialise the unit being assigned to, and the approximate number of blocks
      // required on each unit.
      proc_t currentUnit = 0;

      site_t blocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (mTopologySize));

      // Create an array to monitor whether each block has been assigned yet.
      bool *blockAssigned = new bool[blockCount];

      for (site_t ii = 0; ii < blockCount; ++ii)
      {
        blockAssigned[ii] = false;
      }

      std::vector<BlockLocation> *lCurrentEdge = new std::vector<BlockLocation>;
      std::vector<BlockLocation> *lExpandedEdge = new std::vector<BlockLocation>;

      int lBlockNumber = -1;

      // Domain Decomposition.  Pick a site. Set it to the rank we are
      // looking at. Find its neighbours and put those on the same
      // rank, then find the next-nearest neighbours, etc. until we
      // have a completely joined region, or there are enough fluid
      // sites on the rank.  In the former case, start again at
      // another site. In the latter case, move on to the next rank.
      // Do this until all sites are assigned to a rank. There is a
      // high chance of of all sites on a rank being joined.

      site_t lBlocksOnCurrentProc = 0;

      // Iterate over all blocks.
      for (site_t lBlockCoordI = 0; lBlockCoordI < iGlobLatDat->GetXBlockCount(); lBlockCoordI++)
      {
        for (site_t lBlockCoordJ = 0; lBlockCoordJ < iGlobLatDat->GetYBlockCount(); lBlockCoordJ++)
        {
          for (site_t lBlockCoordK = 0; lBlockCoordK < iGlobLatDat->GetZBlockCount(); lBlockCoordK++)
          {
            // Block number is the number of the block we're currently on.
            lBlockNumber++;

            // If the array of proc rank for each site is NULL, we're on an all-solid block.
            // Alternatively, if this block has already been assigned, move on.
            if (fluidSitesPerBlock[lBlockNumber] == 0)
            {
              unitForEachBlock[lBlockNumber] = -1;
              continue;
            }
            else if (blockAssigned[lBlockNumber])
            {
              continue;
            }

            // Assign this block to the current unit.
            blockAssigned[lBlockNumber] = true;
            unitForEachBlock[lBlockNumber] = currentUnit;

            ++lBlocksOnCurrentProc;

            // Record the location of this initial site.
            lCurrentEdge->clear();
            BlockLocation lNew;
            lNew.i = lBlockCoordI;
            lNew.j = lBlockCoordJ;
            lNew.k = lBlockCoordK;
            lCurrentEdge->push_back(lNew);

            // The subdomain can grow.
            bool lIsRegionGrowing = true;

            // While the region can grow (i.e. it is not bounded by solids or visited
            // sites), and we need more sites on this particular rank.
            while (lBlocksOnCurrentProc < blocksPerUnit && lIsRegionGrowing)
            {
              lExpandedEdge->clear();

              // Sites added to the edge of the mClusters during the iteration.
              lIsRegionGrowing = Expand(lCurrentEdge,
                                        lExpandedEdge,
                                        iGlobLatDat,
                                        fluidSitesPerBlock,
                                        blockAssigned,
                                        currentUnit,
                                        unitForEachBlock,
                                        lBlocksOnCurrentProc,
                                        blocksPerUnit);

              // When the new layer of edge sites has been found, swap the buffers for
              // the current and new layers of edge sites.
              std::vector<BlockLocation> *tempP = lCurrentEdge;
              lCurrentEdge = lExpandedEdge;
              lExpandedEdge = tempP;
            }

            // If we have enough sites, we have finished.
            if (lBlocksOnCurrentProc >= blocksPerUnit)
            {
              blocksOnEachUnit[currentUnit] = lBlocksOnCurrentProc;

              ++currentUnit;

              unassignedBlocks -= lBlocksOnCurrentProc;
              blocksPerUnit = (int) ceil((double) unassignedBlocks / (double) (unitCount
                  - currentUnit));

              lBlocksOnCurrentProc = 0;
            }
            // If not, we have to start growing a different region for the same rank:
            // region expansions could get trapped.

          } // Block co-ord k
        } // Block co-ord j
      } // Block co-ord i

      delete lCurrentEdge;
      delete lExpandedEdge;

      delete[] blockAssigned;
    }

    /**
     * Returns true if the region was expanded.
     *
     * @param edgeBlocks
     * @param iGlobLatDat
     * @param fluidSitesPerBlock
     * @param blockAssigned
     * @param currentUnit
     * @param unitForEachBlock
     * @param blocksOnCurrentProc
     * @param blocksPerUnit
     * @return
     */
    bool LatticeData::GeometryReader::Expand(std::vector<BlockLocation>* edgeBlocks,
                                             std::vector<BlockLocation>* expansionBlocks,
                                             const GlobalLatticeData* iGlobLatDat,
                                             const site_t* fluidSitesPerBlock,
                                             bool* blockAssigned,
                                             proc_t currentUnit,
                                             proc_t* unitForEachBlock,
                                             site_t &blocksOnCurrentUnit,
                                             site_t blocksPerUnit)
    {
      bool lRet = false;

      // For sites on the edge of the domain (sites_a), deal with the neighbours.
      for (unsigned int index_a = 0; index_a < edgeBlocks->size() && blocksOnCurrentUnit
          < blocksPerUnit; index_a++)
      {
        BlockLocation* lNew = &edgeBlocks->operator [](index_a);

        for (unsigned int l = 1; l < D3Q15::NUMVECTORS && blocksOnCurrentUnit < blocksPerUnit; l++)
        {
          // Record neighbour location.
          site_t neigh_i = lNew->i + D3Q15::CX[l];
          site_t neigh_j = lNew->j + D3Q15::CY[l];
          site_t neigh_k = lNew->k + D3Q15::CZ[l];

          // Move on if neighbour is outside the bounding box.
          if (neigh_i == -1 || neigh_i == iGlobLatDat->GetXBlockCount())
            continue;
          if (neigh_j == -1 || neigh_j == iGlobLatDat->GetYBlockCount())
            continue;
          if (neigh_k == -1 || neigh_k == iGlobLatDat->GetZBlockCount())
            continue;

          // Move on if the neighbour is in a block of solids (in which case
          // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid or has already
          // been assigned to a rank (in which case ProcessorRankForEachBlockSite != -1).  ProcessorRankForEachBlockSite
          // was initialized in lbmReadConfig in io.cc.

          site_t neighBlockId = iGlobLatDat->GetBlockIdFromBlockCoords(neigh_i, neigh_j, neigh_k);

          // Don't use this block if it has no fluid sites, or if it has already been assigned to a processor.
          if (fluidSitesPerBlock[neighBlockId] == 0 || blockAssigned[neighBlockId])
          {
            continue;
          }

          // Set the rank for a neighbour and update the fluid site counters.
          blockAssigned[neighBlockId] = true;
          unitForEachBlock[neighBlockId] = currentUnit;
          ++blocksOnCurrentUnit;

          // Neighbour was found, so the region can grow.
          lRet = true;

          // Record the location of the neighbour.
          BlockLocation lNewB;
          lNewB.i = neigh_i;
          lNewB.j = neigh_j;
          lNewB.k = neigh_k;
          expansionBlocks->push_back(lNewB);
        }
      }

      return lRet;
    }

    void LatticeData::GeometryReader::OptimiseDomainDecomposition(const site_t* sitesPerBlock,
                                                                  const unsigned int* bytesPerBlock,
                                                                  const proc_t* procForEachBlock,
                                                                  MPI_File iFile,
                                                                  GlobalLatticeData* bGlobLatDat)
    {
      idxtype* vertexDistribution = new idxtype[mTopologySize + 1];

      /*
       *  Get an array of the site count on each processor.
       */
      {
        site_t* sitesPerProc = new site_t[mTopologySize];
        for (unsigned int ii = 0; ii < mTopologySize; ++ii)
        {
          sitesPerProc[ii] = 0;
        }

        for (unsigned int ii = 0; ii < bGlobLatDat->GetBlockCount(); ++ii)
        {
          if (procForEachBlock[ii] >= 0)
          {
            sitesPerProc[procForEachBlock[ii]] += sitesPerBlock[ii];
          }
        }

        /*
         *  Create vertex distribution array.
         */
        vertexDistribution[0] = 0;

        for (unsigned int ii = 0; ii < mTopologySize; ++ii)
        {
          vertexDistribution[ii + 1] = vertexDistribution[ii] + (int) sitesPerProc[ii];
        }
        delete[] sitesPerProc;
      }

      /*
       *  Create the adjacency count and list for all the local vertices.
       */
      // First, it'll be useful to create an array of the first site index on each block.
      int* firstSiteIndexPerBlock = new int[bGlobLatDat->GetBlockCount()];

      {
        // This needs to have all the blocks with ascending id based on which processor they're on.
        site_t* firstSiteOnProc = new site_t[mTopologySize];
        for (unsigned int ii = 0; ii < mTopologySize; ++ii)
        {
          firstSiteOnProc[ii] = vertexDistribution[ii];
        }

        for (unsigned int ii = 0; ii < bGlobLatDat->GetBlockCount(); ++ii)
        {
          int proc = procForEachBlock[ii];
          if (proc < 0)
          {
            firstSiteIndexPerBlock[ii] = -1;
          }
          else
          {
            firstSiteIndexPerBlock[ii] = (int) firstSiteOnProc[proc];
            firstSiteOnProc[proc] += sitesPerBlock[ii];
          }
        }
        delete[] firstSiteOnProc;
      }

      unsigned int localVertexCount = vertexDistribution[mTopologyRank + 1]
          - vertexDistribution[mTopologyRank];

      idxtype* adjacenciesPerVertex = new idxtype[localVertexCount + 1];
      adjacenciesPerVertex[0] = 0;

      std::vector<idxtype> lAdjacencies;

      unsigned int lFluidVertex = 0;
      int n = -1;
      for (site_t i = 0; i < bGlobLatDat->GetXSiteCount(); i += bGlobLatDat->GetBlockSize())
      {
        for (site_t j = 0; j < bGlobLatDat->GetYSiteCount(); j += bGlobLatDat->GetBlockSize())
        {
          for (site_t k = 0; k < bGlobLatDat->GetZSiteCount(); k += bGlobLatDat->GetBlockSize())
          {
            ++n;

            if (procForEachBlock[n] != mTopologyRank)
            {
              continue;
            }

            BlockData *map_block_p = &bGlobLatDat->Blocks[n];

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            int m = -1;

            // Iterate over sites within the block.
            for (site_t site_i = i; site_i < i + bGlobLatDat->GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + bGlobLatDat->GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + bGlobLatDat->GetBlockSize(); site_k++)
                {
                  ++m;

                  // Continue if it's a solid
                  if (map_block_p->ProcessorRankForEachBlockSite[m] == BIG_NUMBER2)
                  {
                    continue;
                  }

                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Work out positions of neighbours.
                    site_t neigh_i = site_i + D3Q15::CX[l];
                    site_t neigh_j = site_j + D3Q15::CY[l];
                    site_t neigh_k = site_k + D3Q15::CZ[l];

                    if (neigh_i <= 0 || neigh_j <= 0 || neigh_k <= 0
                        || !bGlobLatDat->IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      continue;
                    }

                    // Get the id of the processor which the neighbouring site lies on.
                    const int *proc_id_p = bGlobLatDat->GetProcIdFromGlobalCoords(neigh_i,
                                                                                  neigh_j,
                                                                                  neigh_k);

                    if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
                    {
                      continue;
                    }

                    // We now do some faffery to find out the global fluid site id of this point
                    site_t neighBlockI = neigh_i >> bGlobLatDat->Log2BlockSize;
                    site_t neighBlockJ = neigh_j >> bGlobLatDat->Log2BlockSize;
                    site_t neighBlockK = neigh_k >> bGlobLatDat->Log2BlockSize;

                    site_t neighLocalSiteI = neigh_i - (neighBlockI << bGlobLatDat->Log2BlockSize);
                    site_t neighLocalSiteJ = neigh_j - (neighBlockJ << bGlobLatDat->Log2BlockSize);
                    site_t neighLocalSiteK = neigh_k - (neighBlockK << bGlobLatDat->Log2BlockSize);

                    site_t neighBlockId = bGlobLatDat->GetBlockIdFromBlockCoords(neighBlockI,
                                                                                 neighBlockJ,
                                                                                 neighBlockK);

                    site_t neighGlobalSiteId = firstSiteIndexPerBlock[neighBlockId];

                    site_t localSiteId = ( ( (neighLocalSiteI << bGlobLatDat->Log2BlockSize)
                        + neighLocalSiteJ) << bGlobLatDat->Log2BlockSize) + neighLocalSiteK;

                    for (site_t neighSite = 0; neighSite
                        < bGlobLatDat->GetSitesPerBlockVolumeUnit(); ++neighSite)
                    {
                      if (neighSite == localSiteId)
                      {
                        break;
                      }
                      else if (bGlobLatDat->Blocks[neighBlockId].ProcessorRankForEachBlockSite[neighSite]
                          != BIG_NUMBER2)
                      {
                        ++neighGlobalSiteId;
                      }
                    }

                    lAdjacencies.push_back((idxtype) neighGlobalSiteId);
                  }

                  adjacenciesPerVertex[++lFluidVertex] = (idxtype) lAdjacencies.size();
                }
              }
            }
          }
        }
      }

      if (lFluidVertex != localVertexCount)
      {
        std::cerr << "Encountered a different number of vertices on two different parses: "
            << lFluidVertex << " and " << localVertexCount << "\n";
      }

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

      int weightFlag = 2;
      int numberingFlag = 0;
      int edgesCut = 0;
      idxtype* partitionVector = new idxtype[localVertexCount];
      int options[4] = { 0, 0, 0, 0 };

      int myLowest = vertexDistribution[mTopologyRank];
      int myHighest = vertexDistribution[mTopologyRank + 1] - 1;

      for (unsigned int ii = 0; ii < localVertexCount; ++ii)
      {
        partitionVector[ii] = -1;
      }

      int desiredPartitionSize = mTopologySize;

      int noConstraints = 1;
      idxtype* vertexWeight = new idxtype[localVertexCount];
      for (unsigned int ii = 0; ii < localVertexCount; ++ii)
      {
        vertexWeight[ii] = 1;
      }

      float* domainWeights = new float[desiredPartitionSize];

      for (int ii = 0; ii < desiredPartitionSize; ++ii)
      {
        domainWeights[ii] = 1.0F / ((float) desiredPartitionSize);
      }

      float tolerance = 1.005F;

      log::Logger::Log<log::Debug, log::OnePerCore>("Calling ParMetis");
      ParMETIS_V3_PartKway(vertexDistribution,
                           adjacenciesPerVertex,
                           &lAdjacencies[0],
                           vertexWeight,
                           NULL,
                           &weightFlag,
                           &numberingFlag,
                           &noConstraints,
                           &desiredPartitionSize,
                           domainWeights,
                           &tolerance,
                           options,
                           &edgesCut,
                           partitionVector,
                           &mTopologyComm);

      delete[] domainWeights;
      delete[] vertexWeight;
      delete[] adjacenciesPerVertex;
      delete[] vertexDistribution;

      // Right. Let's count how many sites we're going to have to move.
      int moves = 0;
      std::vector<int> moveData;

      for (int ii = 0; ii <= (myHighest - myLowest); ++ii)
      {
        if (partitionVector[ii] != mTopologyRank)
        {
          ++moves;
          moveData.push_back(myLowest + ii);
          moveData.push_back(partitionVector[ii]);
        }
      }

      delete[] partitionVector;

      // Spread this data around, so all processes now how many moves each process is doing.
      int* allMoves = new int[mTopologySize];

      MPI_Allgather(&moves, 1, MPI_INT, allMoves, 1, MPI_INT, mTopologyComm);

      // Count the total moves.
      int totalMoves = 0;

      for (unsigned int ii = 0; ii < mTopologySize; ++ii)
      {
        totalMoves += allMoves[ii];
      }

      // Now share all the lists of moves.
      MPI_Datatype lMoveType;
      MPI_Type_contiguous(2, MPI_INT, &lMoveType);
      MPI_Type_commit(&lMoveType);

      int* movesList = new int[2 * totalMoves];
      int* offsets = new int[mTopologySize];

      offsets[0] = 0;

      for (unsigned int ii = 1; ii < mTopologySize; ++ii)
      {
        offsets[ii] = offsets[ii - 1] + allMoves[ii - 1];
      }

      MPI_Allgatherv(&moveData[0],
                     moves,
                     lMoveType,
                     movesList,
                     allMoves,
                     offsets,
                     lMoveType,
                     mTopologyComm);

      int* newProcForEachBlock = new int[bGlobLatDat->GetBlockCount()];

      for (unsigned int lBlockNumber = 0; lBlockNumber < bGlobLatDat->GetBlockCount(); ++lBlockNumber)
      {
        newProcForEachBlock[lBlockNumber] = procForEachBlock[lBlockNumber];
      }

      // Hackily, set the proc for each block to be this proc whenever we have a proc on that block under the new scheme.
      unsigned int moveIndex = 0;

      for (unsigned int lFromProc = 0; lFromProc < mTopologySize; ++lFromProc)
      {
        for (int lMoveNumber = 0; lMoveNumber < allMoves[lFromProc]; ++lMoveNumber)
        {
          int fluidSite = movesList[2 * moveIndex];
          int toProc = movesList[2 * moveIndex + 1];
          ++moveIndex;

          if (toProc == mTopologyRank)
          {
            unsigned int fluidSiteBlock = 0;

            while ( (procForEachBlock[fluidSiteBlock] < 0)
                || (firstSiteIndexPerBlock[fluidSiteBlock] > fluidSite)
                || ( (firstSiteIndexPerBlock[fluidSiteBlock] + (int) sitesPerBlock[fluidSiteBlock])
                    <= fluidSite))
            {
              fluidSiteBlock++;
            }

            newProcForEachBlock[fluidSiteBlock] = mTopologyRank;
          }
        }
      }

      MPI_Type_free(&lMoveType);

      ReadInLocalBlocks(iFile, bytesPerBlock, newProcForEachBlock, mTopologyRank, bGlobLatDat);

      // Set the proc rank to what it originally was.
      for (unsigned int lBlock = 0; lBlock < bGlobLatDat->GetBlockCount(); ++lBlock)
      {
        int lOriginalProc = procForEachBlock[lBlock];

        if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite != NULL)
        {
          for (unsigned int lSiteIndex = 0; lSiteIndex < bGlobLatDat->GetSitesPerBlockVolumeUnit(); ++lSiteIndex)
          {
            if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                != BIG_NUMBER2)
            {
              // This should be what it originally was...
              // TODO ugly
              bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                  = (mTopologyRank == mGlobalRank)
                    ? lOriginalProc
                    : (lOriginalProc + 1);
            }
          }
        }
      }

      // Then go through and implement the ParMetis partition.
      moveIndex = 0;

      for (unsigned int lFromProc = 0; lFromProc < mTopologySize; ++lFromProc)
      {
        for (int lMoveNumber = 0; lMoveNumber < allMoves[lFromProc]; ++lMoveNumber)
        {
          int fluidSite = movesList[2 * moveIndex];
          int toProc = movesList[2 * moveIndex + 1];
          ++moveIndex;

          unsigned int fluidSiteBlock = 0;

          while ( (procForEachBlock[fluidSiteBlock] < 0) || (firstSiteIndexPerBlock[fluidSiteBlock]
              > fluidSite) || ( (firstSiteIndexPerBlock[fluidSiteBlock]
              + (int) sitesPerBlock[fluidSiteBlock]) <= fluidSite))
          {
            fluidSiteBlock++;
          }

          if (bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite != NULL)
          {
            int fluidSitesToPass = fluidSite - firstSiteIndexPerBlock[fluidSiteBlock];

            unsigned int lSiteIndex = 0;

            while (true)
            {
              if (bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                  != BIG_NUMBER2)
              {
                fluidSitesToPass--;
              }
              if (fluidSitesToPass < 0)
              {
                break;
              }
              lSiteIndex++;
            }

            // TODO Ugly.
            bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                = (mGlobalRank == mTopologyRank)
                  ? toProc
                  : (toProc + 1);
          }
        }
      }

      delete[] movesList;
      delete[] offsets;
      delete[] firstSiteIndexPerBlock;
      delete[] newProcForEachBlock;
      delete[] allMoves;
    }

    // The header section of the config file contains two adjacent unsigned ints for each block.
    // The first is the number of sites in that block, the second is the number of bytes used by
    // the block in the data file.
    site_t LatticeData::GeometryReader::GetHeaderLength(site_t blockCount) const
    {
      return 2 * 4 * blockCount;
    }

  }
}
