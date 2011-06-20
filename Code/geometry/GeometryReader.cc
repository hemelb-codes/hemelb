#include <math.h>
#include <list>

#include "debug/Debugger.h"
#include "io/XdrMemReader.h"
#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"
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
      // Get the group of all procs.
      MPI_Group lWorldGroup;
      MPI_Comm_group(MPI_COMM_WORLD, &lWorldGroup);

      mParticipateInTopology = !reserveSteeringCore
          || topology::NetworkTopology::Instance()->GetLocalRank() != 0;

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

      int lError;

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      MPI_Info_create(&fileInfo);

      // Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
      std::string accessStyle = "access_style";
      std::string accessStyleValue = "sequential";
      std::string buffering = "collective_buffering";
      std::string bufferingValue = "true";

      MPI_Info_set(fileInfo, &accessStyle[0], &accessStyleValue[0]);
      MPI_Info_set(fileInfo, &buffering[0], &bufferingValue[0]);

      // Open the file.
      lError = MPI_File_open(MPI_COMM_WORLD,
                             &bSimConfig->DataFilePath[0],
                             MPI_MODE_RDONLY,
                             fileInfo,
                             &file);

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
      MPI_File_set_view(file, 0, MPI_CHAR, MPI_CHAR, &lMode[0], fileInfo);

      // Read the file preamble.
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading file preamble");
      ReadPreamble(bLbmParams, bGlobLatDat);

      // Read the file header.
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading file header");

      site_t* sitesPerBlock = new site_t[bGlobLatDat->GetBlockCount()];
      unsigned int* bytesPerBlock = new unsigned int[bGlobLatDat->GetBlockCount()];

      ReadHeader(bGlobLatDat->GetBlockCount(), sitesPerBlock, bytesPerBlock);

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

      // Close the file - only the ranks participating in the topology need to read it again.
      MPI_File_close(&file);

      if (mParticipateInTopology)
      {
        // Reopen in the file just between the nodes in the topology decomposition. Read in blocks
        // local to this node.
        MPI_File_open(mTopologyComm, &bSimConfig->DataFilePath[0], MPI_MODE_RDONLY, fileInfo, &file);

        ReadInLocalBlocks(bGlobLatDat,
                          sitesPerBlock,
                          bytesPerBlock,
                          procForEachBlock,
                          mTopologyRank);
      }

      double lMiddle = util::myClock();

      // Optimise.
      if (mParticipateInTopology)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Beginning domain decomposition optimisation");
        OptimiseDomainDecomposition(sitesPerBlock, bytesPerBlock, procForEachBlock, bGlobLatDat);
        log::Logger::Log<log::Debug, log::OnePerCore>("Ending domain decomposition optimisation");

        MPI_File_close(&file);
      }

      // Finish up - close the file, set the timings, deallocate memory.
      MPI_Info_free(&fileInfo);

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
    void LatticeData::GeometryReader::ReadPreamble(hemelb::lb::LbmParameters * bParams,
                                                   GlobalLatticeData* bGlobalLatticeData)
    {
      // Read in the file preamble into a buffer.
      char lPreambleBuffer[PreambleBytes];

      MPI_Status lStatus;

      MPI_File_read_all(file,
                        lPreambleBuffer,
                        PreambleBytes,
                        MpiDataType(lPreambleBuffer[0]),
                        &lStatus);

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
    void LatticeData::GeometryReader::ReadHeader(site_t iBlockCount,
                                                 site_t* sitesInEachBlock,
                                                 unsigned int* bytesUsedByBlockInDataFile)
    {
      site_t headerByteCount = GetHeaderLength(iBlockCount);
      // Allocate a buffer to read into, then do the reading.
      char* lHeaderBuffer = new char[headerByteCount];

      MPI_Status lStatus;

      MPI_File_read_all(file,
                        lHeaderBuffer,
                        (int) headerByteCount,
                        MpiDataType(lHeaderBuffer[0]),
                        &lStatus);

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
    void LatticeData::GeometryReader::ReadInLocalBlocks(GlobalLatticeData* iGlobLatDat,
                                                        const site_t* sitesPerBlock,
                                                        const unsigned int* bytesPerBlock,
                                                        const proc_t* unitForEachBlock,
                                                        const proc_t localRank)
    {
      // Create a list of which blocks to read in.
      bool* readBlock = new bool[iGlobLatDat->GetBlockCount()];

      DecideWhichBlocksToRead(readBlock, unitForEachBlock, localRank, iGlobLatDat);

      // Set the view and read in.
      MPI_Offset lOffset = PreambleBytes + GetHeaderLength(iGlobLatDat->GetBlockCount());

      // Allocate a buffer to read into.
      const site_t maxBytesPerBlock = (iGlobLatDat->GetSitesPerBlockVolumeUnit()) * (4 * 1 + 8 * (4
          + D3Q15::NUMVECTORS - 1));
      const unsigned int MaxBlocksToReadInOneGo = 100;

      char* buffer = new char[maxBytesPerBlock * MaxBlocksToReadInOneGo];

      // Track the next block we should look at.
      site_t nextBlockToRead = 0;

      // While the next block is valid...
      while (nextBlockToRead < iGlobLatDat->GetBlockCount())
      {
        // ... track the blocks we're going to read, and how many bytes we'll require.
        std::vector<site_t> thisReadBlocks;
        int length = 0;

        // We only want to read valid blocks, up to MaxBlocks blocks, and we only want to read-in
        // consecutive blocks.
        for (; nextBlockToRead < iGlobLatDat->GetBlockCount(); ++nextBlockToRead)
        {
          if (thisReadBlocks.size() == MaxBlocksToReadInOneGo)
          {
            break;
          }

          // If this is a block we're not going to read in...
          if (bytesPerBlock[nextBlockToRead] == 0 || !readBlock[nextBlockToRead])
          {
            // ... and we've not read any blocks yet, add to the start position of the read
            // and keep going until we find a block we do want to read.
            if (thisReadBlocks.size() == 0)
            {
              if (bytesPerBlock[nextBlockToRead] > 0)
              {
                lOffset += bytesPerBlock[nextBlockToRead];
              }

              continue;
            }
            // Alternatively, if we've read some blocks, break because we'd read a
            // non-contiguous set if we continued.

            else
            {
              break;
            }
          }

          // We're going to read this block. Add it to the list, add its length to the total
          // read length.
          thisReadBlocks.push_back(nextBlockToRead);
          length += bytesPerBlock[nextBlockToRead];
        }

        // If the length is going to overflow the buffer print an error.
        if (length > (MaxBlocksToReadInOneGo * maxBytesPerBlock))
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("Read in %i bytes when the longest read was presumed to be: ",
                                                        length,
                                                         (MaxBlocksToReadInOneGo * maxBytesPerBlock));
        }

        // Read the data.
        int error = MPI_File_read_at(file, lOffset, buffer, length, MPI_CHAR, MPI_STATUS_IGNORE);

        // Make sure we don't re-read a section of the file next time.
        lOffset += length;

        // Check for any errors in the file-reading.
        if (error != 0)
        {
          log::Logger::Log<log::Info, log::OnePerCore>("Unable to read file (at section just before block %i)",
                                                       nextBlockToRead);
        }

        // Create an Xdr interpreter.
        io::XdrMemReader lReader(buffer, length);

        unsigned int lCurrentPosition = lReader.GetPosition();

        // Go through the blocks.
        for (std::vector<site_t>::const_iterator it = thisReadBlocks.begin(); it
            < thisReadBlocks.end(); ++it)
        {
          // Let the GlobLatDat read the block.
          iGlobLatDat->ReadBlock(*it, &lReader);

          // If debugging, check we've read the right length.
          if (log::Logger::ShouldDisplay<log::Debug>())
          {
            unsigned int lNewPosition = lReader.GetPosition();

            // Check we've moved the right amount through the file.
            if (lNewPosition != (lCurrentPosition + bytesPerBlock[*it]))
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Was expecting block %i to take up %i bytes but it actually took up %u bytes",
                                                            *it,
                                                             (lNewPosition - lCurrentPosition),
                                                            bytesPerBlock[*it]);
            }

            // Check the bytes for each block fit our expectations.
            if (bytesPerBlock[*it] > maxBytesPerBlock)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Block %i takes up %i bytes but we thought the maximum was %u bytes",
                                                            *it,
                                                            maxBytesPerBlock);
            }

            lCurrentPosition = lNewPosition;
          }
        }
      }

      // If debug-level logging, check that we've read in as many sites as anticipated.
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        for (site_t block = 0; block < iGlobLatDat->GetBlockCount(); ++block)
        {
          // For each block to be read.
          if (bytesPerBlock[block] > 0 && readBlock[block])
          {
            // Count the sites read,
            site_t numSitesRead = 0;
            for (site_t site = 0; site < iGlobLatDat->GetSitesPerBlockVolumeUnit(); ++site)
            {
              if (iGlobLatDat->Blocks[block].ProcessorRankForEachBlockSite[site] != BIG_NUMBER2)
              {
                ++numSitesRead;
              }
            }

            // Compare with the sites we expected to read.
            if (numSitesRead != sitesPerBlock[block])
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Was expecting %i fluid sites on block %i but actually read %i",
                                                            sitesPerBlock[block],
                                                            block,
                                                            numSitesRead);
            }
          }
        }
      }

      delete[] buffer;
      delete[] readBlock;
    }

    /**
     * Compile a list of blocks to be read onto this core, including all the ones we perform
     * LB on, and also any of their neighbouring blocks.
     *
     * NOTE: that the skipping-over of blocks without any fluid sites is dealt with by other
     * code.
     *
     * @param readBlock
     * @param unitForEachBlock
     * @param localRank
     * @param iGlobLatDat
     */
    void LatticeData::GeometryReader::DecideWhichBlocksToRead(bool* readBlock,
                                                              const proc_t* unitForEachBlock,
                                                              const proc_t localRank,
                                                              const GlobalLatticeData* iGlobLatDat)
    {
      // Initialise each block to not being read.
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

            if (unitForEachBlock[lBlockId] != localRank)
            {
              continue;
            }

            // Read in all neighbouring blocks.
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

                  readBlock[lNeighId] = true;
                }
              }
            }
          }
        }
      }
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

      site_t lBlockNumber = -1;

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
              blocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (unitCount
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
                                                                  GlobalLatticeData* bGlobLatDat)
    {
      // Get some arrays that ParMetis needs.
      idxtype* vtxDistribn = new idxtype[mTopologySize + 1];

      GetSiteDistributionArray(vtxDistribn,
                               bGlobLatDat->GetBlockCount(),
                               procForEachBlock,
                               sitesPerBlock);

      idxtype* firstSiteIndexPerBlock = new idxtype[bGlobLatDat->GetBlockCount()];

      GetFirstSiteIndexOnEachBlock(firstSiteIndexPerBlock,
                                   bGlobLatDat->GetBlockCount(),
                                   vtxDistribn,
                                   procForEachBlock,
                                   sitesPerBlock);

      idxtype localVertexCount = vtxDistribn[mTopologyRank + 1] - vtxDistribn[mTopologyRank];

      idxtype* adjacenciesPerVertex = new idxtype[localVertexCount + 1];
      std::vector<idxtype> lAdjacencies;

      GetAdjacencyData(adjacenciesPerVertex,
                       lAdjacencies,
                       localVertexCount,
                       procForEachBlock,
                       firstSiteIndexPerBlock,
                       bGlobLatDat);

      // Call parmetis.
      idxtype* partitionVector = new idxtype[localVertexCount];

      CallParmetis(partitionVector,
                   localVertexCount,
                   vtxDistribn,
                   adjacenciesPerVertex,
                   &lAdjacencies[0]);

      delete[] adjacenciesPerVertex;

      // Convert the ParMetis results into a nice format.
      idxtype* allMoves = new idxtype[mTopologySize];

      idxtype* movesList = GetMovesList(allMoves,
                                        firstSiteIndexPerBlock,
                                        procForEachBlock,
                                        sitesPerBlock,
                                        vtxDistribn,
                                        partitionVector,
                                        bGlobLatDat);

      delete[] firstSiteIndexPerBlock;
      delete[] vtxDistribn;
      delete[] partitionVector;

      // Reread the blocks based on the ParMetis decomposition.
      RereadBlocks(bGlobLatDat, allMoves, movesList, sitesPerBlock, bytesPerBlock, procForEachBlock);

      // Implement the decomposition now that we have read the necessary data.
      ImplementMoves(bGlobLatDat, procForEachBlock, allMoves, movesList);

      delete[] movesList;
      delete[] allMoves;
    }

    // The header section of the config file contains two adjacent unsigned ints for each block.
    // The first is the number of sites in that block, the second is the number of bytes used by
    // the block in the data file.
    site_t LatticeData::GeometryReader::GetHeaderLength(site_t blockCount) const
    {
      return 2 * 4 * blockCount;
    }

    /**
     * Get the cumulative count of sites on each processor.
     *
     * @param vertexDistribn
     * @param blockCount
     * @param procForEachBlock
     * @param sitesPerBlock
     */
    void LatticeData::GeometryReader::GetSiteDistributionArray(idxtype* vertexDistribn,
                                                               const site_t blockCount,
                                                               const proc_t* procForEachBlock,
                                                               const site_t* sitesPerBlock) const
    {
      // Firstly, count the sites per processor.
      for (unsigned int ii = 0; ii < (mTopologySize + 1); ++ii)
      {
        vertexDistribn[ii] = 0;
      }

      for (site_t ii = 0; ii < blockCount; ++ii)
      {
        if (procForEachBlock[ii] >= 0)
        {
          vertexDistribn[1 + procForEachBlock[ii]] += (idxtype) sitesPerBlock[ii];
        }
      }

      // Now make the count cumulative.
      for (unsigned int ii = 0; ii < mTopologySize; ++ii)
      {
        vertexDistribn[ii + 1] += vertexDistribn[ii];
      }
    }

    void LatticeData::GeometryReader::GetFirstSiteIndexOnEachBlock(idxtype* firstSiteIndexPerBlock,
                                                                   const site_t blockCount,
                                                                   const idxtype* vertexDistribution,
                                                                   const proc_t* procForEachBlock,
                                                                   const site_t* sitesPerBlock) const
    {
      // First calculate the lowest site index on each proc - relatively easy.
      site_t* firstSiteOnProc = new site_t[mTopologySize];
      for (unsigned int ii = 0; ii < mTopologySize; ++ii)
      {
        firstSiteOnProc[ii] = vertexDistribution[ii];
      }

      // Now for each block (in ascending order), the smallest site index is the smallest site
      // index on its processor, incremented by the number of sites observed from that processor
      // so far.
      for (site_t ii = 0; ii < blockCount; ++ii)
      {
        proc_t proc = procForEachBlock[ii];
        if (proc < 0)
        {
          firstSiteIndexPerBlock[ii] = -1;
        }
        else
        {
          firstSiteIndexPerBlock[ii] = (idxtype) firstSiteOnProc[proc];
          firstSiteOnProc[proc] += sitesPerBlock[ii];
        }
      }

      // Clean up.
      delete[] firstSiteOnProc;
    }

    void LatticeData::GeometryReader::GetAdjacencyData(idxtype* adjacenciesPerVertex,
                                                       std::vector<idxtype> &adjacencies,
                                                       const idxtype localVertexCount,
                                                       const proc_t* procForEachBlock,
                                                       const idxtype* firstSiteIndexPerBlock,
                                                       const GlobalLatticeData* bGlobLatDat) const
    {
      adjacenciesPerVertex[0] = 0;
      idxtype lFluidVertex = 0;
      site_t n = -1;

      // For each block (counting up by lowest site id)...
      for (site_t i = 0; i < bGlobLatDat->GetXSiteCount(); i += bGlobLatDat->GetBlockSize())
      {
        for (site_t j = 0; j < bGlobLatDat->GetYSiteCount(); j += bGlobLatDat->GetBlockSize())
        {
          for (site_t k = 0; k < bGlobLatDat->GetZSiteCount(); k += bGlobLatDat->GetBlockSize())
          {
            ++n;

            // ... considering only the ones which live on this proc...
            if (procForEachBlock[n] != mTopologyRank)
            {
              continue;
            }

            BlockData *map_block_p = &bGlobLatDat->Blocks[n];

            // ... and only those with fluid sites...
            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            site_t m = -1;

            // ... iterate over sites within the block...
            for (site_t site_i = i; site_i < i + bGlobLatDat->GetBlockSize(); site_i++)
            {
              for (site_t site_j = j; site_j < j + bGlobLatDat->GetBlockSize(); site_j++)
              {
                for (site_t site_k = k; site_k < k + bGlobLatDat->GetBlockSize(); site_k++)
                {
                  ++m;

                  // ... only looking at non-solid sites...
                  if (map_block_p->ProcessorRankForEachBlockSite[m] == BIG_NUMBER2)
                  {
                    continue;
                  }

                  // ... for each lattice direction...
                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // ... which leads to a valid neighbouring site...
                    site_t neigh_i = site_i + D3Q15::CX[l];
                    site_t neigh_j = site_j + D3Q15::CY[l];
                    site_t neigh_k = site_k + D3Q15::CZ[l];

                    if (neigh_i <= 0 || neigh_j <= 0 || neigh_k <= 0
                        || !bGlobLatDat->IsValidLatticeSite(neigh_i, neigh_j, neigh_k))
                    {
                      continue;
                    }

                    // ... (that is actually being simulated and not a solid)...
                    const proc_t* proc_id_p = bGlobLatDat->GetProcIdFromGlobalCoords(neigh_i,
                                                                                     neigh_j,
                                                                                     neigh_k);

                    if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
                    {
                      continue;
                    }

                    // ... get some info about the position of the neighbouring site.
                    site_t neighBlockI = neigh_i >> bGlobLatDat->Log2BlockSize;
                    site_t neighBlockJ = neigh_j >> bGlobLatDat->Log2BlockSize;
                    site_t neighBlockK = neigh_k >> bGlobLatDat->Log2BlockSize;

                    site_t neighLocalSiteI = neigh_i - (neighBlockI << bGlobLatDat->Log2BlockSize);
                    site_t neighLocalSiteJ = neigh_j - (neighBlockJ << bGlobLatDat->Log2BlockSize);
                    site_t neighLocalSiteK = neigh_k - (neighBlockK << bGlobLatDat->Log2BlockSize);

                    site_t neighBlockId = bGlobLatDat->GetBlockIdFromBlockCoords(neighBlockI,
                                                                                 neighBlockJ,
                                                                                 neighBlockK);

                    // Now get the local id of the neighbour on its block,
                    site_t localSiteId = ( ( (neighLocalSiteI << bGlobLatDat->Log2BlockSize)
                        + neighLocalSiteJ) << bGlobLatDat->Log2BlockSize) + neighLocalSiteK;

                    // calculate the site's id over the whole geometry,
                    site_t neighGlobalSiteId = firstSiteIndexPerBlock[neighBlockId];

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

                    // then add this to the list of adjacencies.
                    adjacencies.push_back((idxtype) neighGlobalSiteId);
                  }

                  // The cumulative count of adjacencies for this vertex is equal to the total
                  // number of adjacencies we've entered.
                  // NOTE: The prefix operator is correct here because
                  // the array has a leading 0 not relating to any site.
                  adjacenciesPerVertex[++lFluidVertex] = (idxtype) adjacencies.size();
                }
              }
            }
          }
        }
      }

      // Perform a debugging test if running at the appropriate log level.
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        if (lFluidVertex != localVertexCount)
        {
          std::cerr << "Encountered a different number of vertices on two different parses: "
              << lFluidVertex << " and " << localVertexCount << "\n";
        }
      }
    }

    void LatticeData::GeometryReader::CallParmetis(idxtype* partitionVector,
                                                   idxtype localVertexCount,
                                                   idxtype* vtxDistribn,
                                                   idxtype* adjacenciesPerVertex,
                                                   idxtype* adjacencies)
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
      for (idxtype ii = 0; ii < localVertexCount; ++ii)
      {
        partitionVector[ii] = mTopologyRank;
      }

      // Weight all vertices evenly.
      idxtype* vertexWeight = new idxtype[localVertexCount];
      for (idxtype ii = 0; ii < localVertexCount; ++ii)
      {
        vertexWeight[ii] = 1;
      }

      // Set the weights of each partition to be even, and to sum to 1.
      idxtype desiredPartitionSize = mTopologySize;

      float* domainWeights = new float[desiredPartitionSize];
      for (idxtype ii = 0; ii < desiredPartitionSize; ++ii)
      {
        domainWeights[ii] = 1.0F / ((float) desiredPartitionSize);
      }

      // A bunch of values ParMetis needs.
      idxtype noConstraints = 1;
      idxtype weightFlag = 2;
      idxtype numberingFlag = 0;
      idxtype edgesCut = 0;
      idxtype options[4] = { 0, 0, 0, 0 };

      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        // Specify that some options are set and that we should
        // debug everything.
        options[0] = 1;
        options[1] = 1 | 2 | 4 | 8 | 32 | 64;
      }

      float tolerance = 1.001F;

      // If we're using debugging logs, check that the arguments are consistent across all cores.
      // To verify: vtxDistribn, adjacenciesPerVertex, adjacencies
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        // vtxDistribn should be the same on all cores.
        idxtype* vtxDistribnRecv = new idxtype[mTopologySize + 1];

        MPI_Allreduce(vtxDistribn,
                      vtxDistribnRecv,
                      mTopologySize + 1,
                      MpiDataType(vtxDistribnRecv[0]),
                      MPI_MIN,
                      mTopologyComm);

        for (unsigned int ii = 0; ii < mTopologySize + 1; ++ii)
        {
          if (vtxDistribn[ii] != vtxDistribnRecv[ii])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("vtxDistribn[%i] was %li but at least one other core had it as %li.",
                                                          ii,
                                                          vtxDistribn[ii],
                                                          vtxDistribnRecv[ii]);
          }
        }

        delete[] vtxDistribnRecv;

        // Create an array of lists to store all of this node's adjacencies, arranged by the
        // proc the adjacent vertex is on.
        std::list<idxtype>* adjByNeighProc = new std::list<idxtype>[mTopologySize];
        for (proc_t ii = 0; ii < (proc_t) mTopologySize; ++ii)
        {
          adjByNeighProc[ii] = std::list<idxtype>();
        }

        // The adjacency data should correspond across all cores.
        // This is likely to be VERY time-consuming.
        for (idxtype index = 0; index < localVertexCount; ++index)
        {
          idxtype vertex = vtxDistribn[mTopologyRank] + index;

          // Iterate over each adjacency (of each vertex).
          for (idxtype adjNumber = 0; adjNumber < (adjacenciesPerVertex[index + 1]
              - adjacenciesPerVertex[index]); ++adjNumber)
          {
            idxtype adjacentVertex = adjacencies[adjacenciesPerVertex[index] + adjNumber];
            proc_t adjacentProc = -1;

            // Calculate the proc of the neighbouring vertex.
            for (proc_t ii = 0; ii < (proc_t) mTopologySize; ++ii)
            {
              if (vtxDistribn[ii] <= adjacentVertex && vtxDistribn[ii + 1] > adjacentVertex)
              {
                adjacentProc = ii;
                break;
              }
            }

            // If it doesn't appear to belong to any proc, something's wrong.
            if (adjacentProc == -1)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("The vertex %li has a neighbour %li which doesn\'t appear to live on any processor.",
                                                            vertex,
                                                            adjacentVertex);
              continue;
            }

            // Store the data if it does belong to a proc.
            adjByNeighProc[adjacentProc].push_back(adjacentVertex);
            adjByNeighProc[adjacentProc].push_back(vertex);
          }
        }

        // Now spread and compare the adjacency information. Smaller ranks send data to larger
        // ranks which receive the data and compare it.
        for (proc_t neigh = 0; neigh < (proc_t) mTopologySize; ++neigh)
        {
          if (neigh > mTopologyRank)
          {
            // Send the array length.
            idxtype count = adjByNeighProc[neigh].size();
            MPI_Send(&count, 1, MpiDataType(count), neigh, 42, mTopologyComm);

            // Create a sendable array (std::lists aren't organised in a sendable format).
            idxtype* sendArray = new idxtype[count];

            unsigned int ii = 0;

            for (std::list<idxtype>::iterator it = adjByNeighProc[neigh].begin(); it
                != adjByNeighProc[neigh].end(); ++it)

            {
              sendArray[ii] = *it;
              ++ii;
            }

            // Send the data to the neighbour.
            MPI_Send(sendArray, (int) count, MpiDataType<idxtype> (), neigh, 43, mTopologyComm);

            delete[] sendArray;

            // Sending arrays don't perform comparison.
            continue;
          }

          idxtype recvCount;
          idxtype* recvData;

          // If this is a greater rank number than the neighbour, receive the data.
          if (neigh < mTopologyRank)
          {
            MPI_Recv(&recvCount,
                     1,
                     MpiDataType(recvCount),
                     neigh,
                     42,
                     mTopologyComm,
                     MPI_STATUS_IGNORE);

            recvData = new idxtype[recvCount];

            MPI_Recv(recvData,
                     (int) recvCount,
                     MpiDataType<idxtype> (),
                     neigh,
                     43,
                     mTopologyComm,
                     MPI_STATUS_IGNORE);
          }
          // Neigh == mTopologyRank, i.e. neighbouring vertices on the same proc
          // Duplicate the data.
          else
          {
            recvCount = adjByNeighProc[neigh].size();
            recvData = new idxtype[recvCount];

            int ii = 0;
            for (std::list<idxtype>::iterator it = adjByNeighProc[neigh].begin(); it
                != adjByNeighProc[neigh].end(); ++it)

            {
              recvData[ii] = *it;
              ++ii;
            }

          }

          // Now we compare. First go through the received data which is ordered as (adjacent
          // vertex, vertex) wrt the neighbouring proc.
          for (idxtype ii = 0; ii < recvCount; ii += 2)
          {
            bool found = false;
            std::list<idxtype>::iterator it = adjByNeighProc[neigh].begin();

            // Go through each neighbour we know about on this proc, and check whether it
            // matches the current received neighbour-data.
            while (it != adjByNeighProc[neigh].end())
            {
              std::list<idxtype>::iterator initial = it;

              idxtype recvAdj = *it;
              it++;
              idxtype recvAdj2 = *it;
              it++;

              if (recvData[ii] == recvAdj2 && recvData[ii + 1] == recvAdj)
              {
                adjByNeighProc[neigh].erase(initial, it);
                found = true;
                break;
              }
            }

            // No neighbour data on this proc matched the data received.
            if (!found)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Neighbour proc %i had adjacency (%li,%li) that wasn't present on this processor.",
                                                            neigh,
                                                            recvData[ii],
                                                            recvData[ii + 1]);
            }
          }

          // The local store of adjacencies should now be empty, if there was complete matching.
          std::list<idxtype>::iterator it = adjByNeighProc[neigh].begin();
          while (it != adjByNeighProc[neigh].end())
          {
            idxtype adj1 = *it;
            ++it;
            idxtype adj2 = *it;
            ++it;

            // We had neighbour-data on this proc that didn't match that received.
            log::Logger::Log<log::Debug, log::OnePerCore>("The local processor has adjacency (%li,%li) that isn't present on neighbouring processor %i.",
                                                          adj1,
                                                          adj2,
                                                          neigh);

          }

          delete[] recvData;
        }

        delete[] adjByNeighProc;
      }

      log::Logger::Log<log::Debug, log::OnePerCore>("Calling ParMetis");
      ParMETIS_V3_PartKway(vtxDistribn,
                           adjacenciesPerVertex,
                           adjacencies,
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
      log::Logger::Log<log::Debug, log::OnePerCore>("ParMetis returned.");

      delete[] domainWeights;
      delete[] vertexWeight;
    }

    /**
     * Returns a list of the fluid sites to be moved.
     *
     * NOTE: This function's return value is a dynamically-allocated array of all the moves to be
     * performed, ordered by (origin processor [with a count described by the content of the first
     * parameter], site id on the origin processor). The contents of the array are contiguous
     * triplets of ints: (block id, site id on block, destination rank).
     *
     * @param movesFromEachProc
     * @param vtxDistribn
     * @param partitionVector
     * @return
     */
    idxtype* LatticeData::GeometryReader::GetMovesList(idxtype* movesFromEachProc,
                                                       const idxtype* firstSiteIndexPerBlock,
                                                       const proc_t* procForEachBlock,
                                                       const site_t* sitesPerBlock,
                                                       const idxtype* vtxDistribn,
                                                       const idxtype* partitionVector,
                                                       const GlobalLatticeData* bGlobLatDat)
    {
      // Right. Let's count how many sites we're going to have to move. Count the local number of
      // sites to be moved, and collect the site id and the destination processor.
      std::vector<idxtype> moveData;

      const idxtype myLowest = vtxDistribn[mTopologyRank];
      const idxtype myHighest = vtxDistribn[mTopologyRank + 1] - 1;

      // For each local fluid site...
      for (idxtype ii = 0; ii <= (myHighest - myLowest); ++ii)
      {
        // ... if it's going elsewhere...
        if (partitionVector[ii] != mTopologyRank)
        {
          // ... get it's id on the local processor...
          idxtype localFluidSiteId = myLowest + ii;

          // ... find out which block it's on, by going through all blocks until we find one
          // with firstIdOnBlock <= localFluidSiteId < (firstIdOnBlock + sitesOnBlock)...
          idxtype fluidSiteBlock = 0;

          while ( (procForEachBlock[fluidSiteBlock] < 0) || (firstSiteIndexPerBlock[fluidSiteBlock]
              > localFluidSiteId) || ( (firstSiteIndexPerBlock[fluidSiteBlock]
              + sitesPerBlock[fluidSiteBlock]) <= localFluidSiteId))
          {
            fluidSiteBlock++;
          }

          // ... and find its site id within that block. Start by working out how many fluid sites
          // we have to pass before we arrive at the fluid site we're after...
          idxtype fluidSitesToPass = localFluidSiteId - firstSiteIndexPerBlock[fluidSiteBlock];
          idxtype lSiteIndex = 0;

          while (true)
          {
            // ... then keep going through the sites on the block until we've passed as many fluid
            // sites as we need to.
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

          // The above code could go wrong, so in debug logging mode, we do some extra tests.
          if (log::Logger::ShouldDisplay<log::Debug>())
          {
            // If we've ended up on an impossible block, or one that doesn't live on this rank,
            // inform the user.
            if (fluidSiteBlock >= bGlobLatDat->GetBlockCount() || procForEachBlock[fluidSiteBlock]
                != mTopologyRank)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Partition element %i wrongly assigned to block %u of %i (block on processor %i)",
                                                            ii,
                                                            fluidSiteBlock,
                                                            bGlobLatDat->GetBlockCount(),
                                                            procForEachBlock[fluidSiteBlock]);
            }

            // Similarly, if we've ended up with an impossible site index, or a solid site,
            // print an error message.
            if (lSiteIndex >= bGlobLatDat->GetSitesPerBlockVolumeUnit()
                || bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                    == BIG_NUMBER2)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Partition element %i wrongly assigned to site %u of %i (block %i%s)",
                                                            ii,
                                                            lSiteIndex,
                                                            sitesPerBlock[fluidSiteBlock],
                                                            fluidSiteBlock,
                                                            bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                                                                == BIG_NUMBER2
                                                              ? " and site is solid"
                                                              : "");
            }
          }

          // Add the block, site and destination rank to our move list.
          moveData.push_back(fluidSiteBlock);
          moveData.push_back(lSiteIndex);
          moveData.push_back(partitionVector[ii]);
        }
      }

      // Spread the move-count data around, so all processes now how many moves each process is
      // doing.
      idxtype moves = moveData.size() / 3;
      MPI_Allgather(&moves,
                    1,
                    MpiDataType(moves),
                    movesFromEachProc,
                    1,
                    MpiDataType(movesFromEachProc[0]),
                    mTopologyComm);

      // Count the total moves.
      idxtype totalMoves = 0;

      for (unsigned int ii = 0; ii < mTopologySize; ++ii)
      {
        totalMoves += movesFromEachProc[ii];
      }

      // Now share all the lists of moves - create a MPI type...
      MPI_Datatype lMoveType;
      MPI_Type_contiguous(3, MpiDataType<idxtype> (), &lMoveType);
      MPI_Type_commit(&lMoveType);

      // ... create a destination array...
      idxtype* movesList = new idxtype[3 * totalMoves];

      // ... create an array of offsets into the destination array for each rank...
      int* offsets = new int[mTopologySize];

      offsets[0] = 0;
      for (unsigned int ii = 1; ii < mTopologySize; ++ii)
      {
        offsets[ii] = (int) (offsets[ii - 1] + movesFromEachProc[ii - 1]);
      }

      {
        int* procMovesInt = new int[mTopologySize];

        for (proc_t ii = 0; ii < (proc_t) mTopologySize; ++ii)
        {
          procMovesInt[ii] = (int) movesFromEachProc[ii];
        }

        // ... use MPI to gather the data...
        MPI_Allgatherv(&moveData[0],
                       (int) moves,
                       lMoveType,
                       movesList,
                       procMovesInt,
                       offsets,
                       lMoveType,
                       mTopologyComm);

        delete[] procMovesInt;

      }

      // ... clean up...
      MPI_Type_free(&lMoveType);
      delete[] offsets;

      // ... and return the list of moves.
      return movesList;
    }

    void LatticeData::GeometryReader::RereadBlocks(GlobalLatticeData* bGlobLatDat,
                                                   const idxtype* movesPerProc,
                                                   const idxtype* movesList,
                                                   const site_t* sitesPerBlock,
                                                   const unsigned int* bytesPerBlock,
                                                   const int* procForEachBlock)
    {
      // Initialise the array (of which proc each block belongs to) to what it was before.
      int* newProcForEachBlock = new int[bGlobLatDat->GetBlockCount()];

      for (site_t lBlockNumber = 0; lBlockNumber < bGlobLatDat->GetBlockCount(); ++lBlockNumber)
      {
        newProcForEachBlock[lBlockNumber] = procForEachBlock[lBlockNumber];
      }

      // Set the proc for each block to be the current proc whenever a site on that block is
      // going to be moved to the current proc.
      idxtype moveIndex = 0;

      for (unsigned int lFromProc = 0; lFromProc < mTopologySize; ++lFromProc)
      {
        for (idxtype lMoveNumber = 0; lMoveNumber < movesPerProc[lFromProc]; ++lMoveNumber)
        {
          idxtype block = movesList[3 * moveIndex];
          idxtype toProc = movesList[3 * moveIndex + 2];
          ++moveIndex;

          if (toProc == (idxtype) mTopologyRank)
          {
            newProcForEachBlock[block] = mTopologyRank;
          }
        }
      }

      // Reread the blocks into the GlobalLatticeData now.
      ReadInLocalBlocks(bGlobLatDat,
                        sitesPerBlock,
                        bytesPerBlock,
                        newProcForEachBlock,
                        mTopologyRank);

      // Clean up.
      delete[] newProcForEachBlock;
    }

    void LatticeData::GeometryReader::ImplementMoves(GlobalLatticeData* bGlobLatDat,
                                                     const proc_t* procForEachBlock,
                                                     const idxtype* movesFromEachProc,
                                                     const idxtype* movesList) const
    {
      // First all, set the proc rank for each site to what it originally was before
      // domain decomposition optimisation. Go through each block...
      for (site_t lBlock = 0; lBlock < bGlobLatDat->GetBlockCount(); ++lBlock)
      {
        // If this proc has owned a fluid site on this block either before or after optimisation,
        // the following will be non-null.
        if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite != NULL)
        {
          // Get the original proc for that block.
          proc_t lOriginalProc = procForEachBlock[lBlock];

          // For each site on that block...
          for (site_t lSiteIndex = 0; lSiteIndex < bGlobLatDat->GetSitesPerBlockVolumeUnit(); ++lSiteIndex)
          {
            // ... if the site is non-solid...
            if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                != BIG_NUMBER2)
            {
              // ... set its rank to be the rank it had before optimisation.
              bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                  = ConvertTopologyRankToGlobalRank(lOriginalProc);
            }
          }
        }
      }

      // Now implement the moves suggested by parmetis.
      idxtype moveIndex = 0;

      // For each source proc, go through as many moves as it had.
      for (unsigned int lFromProc = 0; lFromProc < mTopologySize; ++lFromProc)
      {
        for (idxtype lMoveNumber = 0; lMoveNumber < movesFromEachProc[lFromProc]; ++lMoveNumber)
        {
          // For each move, get the block, site and destination proc.
          idxtype block = movesList[3 * moveIndex];
          idxtype site = movesList[3 * moveIndex + 1];
          idxtype toProc = movesList[3 * moveIndex + 2];

          // Only implement the move if we have read that block's data.
          if (bGlobLatDat->Blocks[block].ProcessorRankForEachBlockSite != NULL)
          {
            // Some logging code - the unmodified rank for each move's site should equal
            // lFromProc.
            if (log::Logger::ShouldDisplay<log::Debug>())
            {
              if (bGlobLatDat->Blocks[block].ProcessorRankForEachBlockSite[site]
                  != ConvertTopologyRankToGlobalRank((proc_t) lFromProc))
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Block %" IDXTYPE_FORMAT ", site %" IDXTYPE_FORMAT " from move %u was originally on proc %i, not proc %u.",
                                                              block,
                                                              site,
                                                              moveIndex,
                                                              bGlobLatDat->Blocks[block].ProcessorRankForEachBlockSite[site],
                                                              lFromProc);
              }
            }

            // Implement the move.
            bGlobLatDat->Blocks[block].ProcessorRankForEachBlockSite[site]
                = ConvertTopologyRankToGlobalRank((proc_t) toProc);
          }

          ++moveIndex;
        }
      }
    }

    proc_t LatticeData::GeometryReader::ConvertTopologyRankToGlobalRank(proc_t topologyRank) const
    {
      // If the global rank is not equal to the topology rank, we are not using rank 0 for
      // LBM.
      return (topology::NetworkTopology::Instance()->GetLocalRank() == mTopologyRank)
        ? topologyRank
        : (topologyRank + 1);
    }

    void LatticeData::GeometryReader::CreateFileReadType(MPI_Datatype* dataType,
                                                         const site_t blockCount,
                                                         const bool* readBlock,
                                                         const unsigned int* bytesPerBlock) const
    {
      // Create vectors for each of the things we'll need to give to MPI_Type_create_struct
      std::vector<MPI_Datatype> baseTypes;
      std::vector<MPI_Aint> displacements;
      std::vector<int> counts;

      int currentDisplacement = 0;

      // For each block, record the type (byte), number of bytes, and distance from the start
      // of the file.
      for (site_t block = 0; block < blockCount; ++block)
      {
        if (readBlock[block])
        {
          baseTypes.push_back(MPI_CHAR);
          counts.push_back(bytesPerBlock[block]);
          displacements.push_back(currentDisplacement);
        }

        currentDisplacement += bytesPerBlock[block];
      }

      // Create the type.
      MPI_Type_create_struct((int) baseTypes.size(),
                             &counts[0],
                             &displacements[0],
                             &baseTypes[0],
                             dataType);
    }

  }
}
