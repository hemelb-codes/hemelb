#include "TopologyReader.h"
#include "io/XdrMemReader.h"
#include "lb/GlobalLatticeData.h"
extern "C"
{
#include "parmetis/parmetislib.h"
}
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace topology
  {

    TopologyReader::TopologyReader()
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &mRank);
      MPI_Comm_size(MPI_COMM_WORLD, &mSize);
      MPI_Comm_group(MPI_COMM_WORLD, &mGroup);

      if (mSize == 1)
      {
        MPI_Group lNewGroup;
        MPI_Comm_dup(MPI_COMM_WORLD, &mCommunicator);
        int lExclusions[1];
        MPI_Group_excl(mGroup, 0, lExclusions, &lNewGroup);
        mGroup = lNewGroup;
      }
      else
      {
        MPI_Group lNewGroup;
        int lExclusions[1] = { 0 };
        MPI_Group_excl(mGroup, 1, lExclusions, &lNewGroup);
        mGroup = lNewGroup;

        MPI_Comm_create(MPI_COMM_WORLD, mGroup, &mCommunicator);
      }
    }

    TopologyReader::~TopologyReader()
    {
      MPI_Group_free(&mGroup);
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
     */
    void TopologyReader::ReadPreamble(MPI_File xiFile,
                                      hemelb::lb::LbmParameters * bParams,
                                      hemelb::lb::GlobalLatticeData* bGlobalLatticeData)
    {
      std::string lMode = "native";

      MPI_File_set_view(xiFile, 0, MPI_BYTE, MPI_BYTE, &lMode[0], MPI_INFO_NULL);

      // The config file starts with:
      // * 1 unsigned int for stress type
      // * 3 unsigned ints for the number of blocks in the x, y, z directions
      // * 1 unsigned int for the block size (number of sites along one edge of a block)
      static const int PreambleBytes = 20;

      char lPreambleBuffer[PreambleBytes];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lPreambleBuffer, PreambleBytes, MPI_BYTE, &lStatus);

      hemelb::io::XdrReader preambleReader = hemelb::io::XdrMemReader(lPreambleBuffer,
                                                                      PreambleBytes);

      // Variables we'll read.
      unsigned int stressType, blocksX, blocksY, blocksZ, blockSize;

      preambleReader.readUnsignedInt(stressType);
      preambleReader.readUnsignedInt(blocksX);
      preambleReader.readUnsignedInt(blocksY);
      preambleReader.readUnsignedInt(blocksZ);
      preambleReader.readUnsignedInt(blockSize);

      bParams->StressType = (lb::StressTypes) stressType;

      bGlobalLatticeData->SetBasicDetails(blocksX, blocksY, blocksZ, blockSize);
    }

    /**
     * Read the header section, with minimal information about each block.
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
    void TopologyReader::ReadHeader(MPI_File xiFile,
                                    unsigned int iBlockCount,
                                    unsigned int* sitesInEachBlock,
                                    unsigned int* bytesUsedByBlockInDataFile)
    {
      // The header section of the config file contains two adjacent unsigned ints for each block.
      // The first is the number of sites in that block, the second is the number of bytes used by
      // the block in the data file.

      unsigned int headerByteCount = 2 * 4 * iBlockCount;

      char* lHeaderBuffer = new char[headerByteCount];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lHeaderBuffer, headerByteCount, MPI_BYTE, &lStatus);

      hemelb::io::XdrReader preambleReader = hemelb::io::XdrMemReader(lHeaderBuffer,
                                                                      headerByteCount);

      for (unsigned int ii = 0; ii < iBlockCount; ii++)
      {
        preambleReader.readUnsignedInt(sitesInEachBlock[ii]);
        preambleReader.readUnsignedInt(bytesUsedByBlockInDataFile[ii]);
      }

      delete[] lHeaderBuffer;
    }

    void TopologyReader::GetInitialSiteDistribution(unsigned long oFirstBlockIdForEachProc[],
                                                    unsigned long oFirstSiteIdForEachProc[],
                                                    unsigned long oNumberSitesPerProc[],
                                                    unsigned long oTotalSiteCount,
                                                    const int iNonSolidSitesPerBlock[],
                                                    const hemelb::lb::GlobalLatticeData & iGlobLatDat)
    {
      oTotalSiteCount = 0;

      for (unsigned int ii = 0; ii < iGlobLatDat.GetBlockCount(); ii++)
      {
        oTotalSiteCount += iNonSolidSitesPerBlock[ii];
      }

      long lApproxPerProc = oTotalSiteCount / mSize;

      int lRoundedProcessors = oTotalSiteCount - (lApproxPerProc * mSize);

      long lCurrentBlockId = 0;
      long lCurrentSiteId = 0;

      for (int ii = 0; ii < mSize; ii++)
      {
        oNumberSitesPerProc[ii] = ( (mSize - lRoundedProcessors) <= ii)
          ? (lApproxPerProc + 1)
          : lApproxPerProc;

        oFirstBlockIdForEachProc[ii] = lCurrentBlockId;
        oFirstSiteIdForEachProc[ii] = lCurrentSiteId;

        int lSitesToPass = oNumberSitesPerProc[ii];

        for (; (lSitesToPass > 0) && (lCurrentBlockId <= iGlobLatDat.GetBlockCount());)
        {
          int lSitesLeftOnCurrentBlock = iNonSolidSitesPerBlock[lCurrentBlockId] - lCurrentSiteId;

          if (lSitesLeftOnCurrentBlock > lSitesToPass)
          {
            lCurrentSiteId += (lSitesLeftOnCurrentBlock - lSitesToPass);
            lSitesToPass = 0;
            break;
          }
          else
          {
            lCurrentSiteId = 0;
            lCurrentBlockId = 0;
            lSitesToPass -= lSitesLeftOnCurrentBlock;
          }
        }
      }
    }

    void TopologyReader::OptimiseDomainDecomposition(const unsigned long iFirstBlockIdForEachProc[],
                                                     const unsigned long iFirstSiteIdForEachProc[],
                                                     const unsigned long iNumberSitesPerProc[],
                                                     const unsigned long iTotalSiteCount,
                                                     const int iNonSolidSitesPerBlock[],
                                                     const hemelb::lb::GlobalLatticeData & iGlobLatDat)
    {
      throw "Not yet implemented";

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

      unsigned long * lCumulativeSitesPerProc = new unsigned long[mSize + 1];
      lCumulativeSitesPerProc[0] = 0;
      for (int ii = 0; ii < mSize; ii++)
      {
        lCumulativeSitesPerProc[ii + 1] = lCumulativeSitesPerProc[ii] + iNumberSitesPerProc[ii];
      }

      //TODO

      //  ParMETIS_RefineKway();
    }

    void TopologyReader::ReadInBlocks(const unsigned long iFirstBlockIdForEachProc[],
                                      const unsigned long iFirstSiteNumberForEachProc[],
                                      const unsigned long iNumberSitesPerProc[],
                                      const unsigned long iTotalSiteCount,
                                      const int iNonSolidSitesPerBlock[],
                                      const hemelb::lb::GlobalLatticeData & iGlobLatDat)
    {
      bool * lReadBlock = new bool[iGlobLatDat.GetBlockCount()];

      for (unsigned int ii = 0; ii < iGlobLatDat.GetBlockCount(); ii++)
      {
        lReadBlock[ii] = false;
      }

      // The largest block we need to read up to is
      // - the last block if this is the last processor
      // - the next proc's first block - 1 if the next block is reading from the very first site on its block.
      // - the next proc's first block if the next block starts reading partway through the site list.
      int lUpperLimit = (mRank == (mSize - 1)
        ? (iGlobLatDat.GetBlockCount() - 1)
        : ( (iFirstSiteNumberForEachProc[mRank + 1] == 0)
          ? iFirstBlockIdForEachProc[mRank + 1] - 1
          : iFirstBlockIdForEachProc[mRank + 1]));

      // For each block that we're computing on, make sure it and the 26 other blocks in the cube around it will
      // be read by this processor.
      for (int ii = iFirstBlockIdForEachProc[mRank]; ii <= lUpperLimit; ii++)
      {
        for (int deltaX = -1; deltaX <= 1; deltaX++)
        {
          for (int deltaY = -1; deltaY <= 1; deltaY++)
          {
            for (int deltaZ = -1; deltaZ <= 1; deltaZ++)
            {
              int lPutativeBlock = ii + ( (deltaX * iGlobLatDat.GetYBlockCount()) + deltaY)
                  * iGlobLatDat.GetZBlockCount() + deltaZ;

              if (lPutativeBlock >= 0 && lPutativeBlock < (int) iGlobLatDat.GetBlockCount())
              {
                lReadBlock[lPutativeBlock] = true;
              }
            }
          }
        }
      }

      // Now we actually do the reading.

    }

    void TopologyReader::LoadAndDecompose(lb::GlobalLatticeData* bGlobalLatticeData,
                                          int &totalFluidSites,
                                          unsigned int siteMins[3],
                                          unsigned int siteMaxes[3],
                                          bool iReserveSteeringCore,
                                          NetworkTopology* bNetTop,
                                          lb::LbmParameters* bLbmParams,
                                          SimConfig* bSimConfig,
                                          double* oReadTime,
                                          double* oDecomposeTime)
    {
      double lStart = MPI_Wtime();

      MPI_File lFile;
      int lError;

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      lError = MPI_File_open(MPI_COMM_WORLD, &bSimConfig->DataFilePath[0], MPI_MODE_RDONLY,
                             MPI_INFO_NULL, &lFile);

      if (lError != 0)
      {
        fprintf(stderr, "Unable to open file %s [rank %i], exiting\n",
                bSimConfig->DataFilePath.c_str(), mRank);
        fflush(0x0);
        exit(0x0);
      }
      else
      {
        fprintf(stderr, "Opened config file %s [rank %i]\n", bSimConfig->DataFilePath.c_str(),
                mRank);
      }
      fflush(NULL);

      ReadPreamble(lFile, bLbmParams, bGlobalLatticeData);

      unsigned int* sitesPerBlock = new unsigned int[bGlobalLatticeData->GetBlockCount()];
      unsigned int* bytesPerBlock = new unsigned int[bGlobalLatticeData->GetBlockCount()];

      ReadHeader(lFile, bGlobalLatticeData->GetBlockCount(), sitesPerBlock, bytesPerBlock);

      ReadAllBlocks(bGlobalLatticeData, bytesPerBlock, totalFluidSites, siteMins, siteMaxes, lFile);
      /*
       * Will eventually become the optimised version...

       unsigned long * lFirstBlockIdForEachProc = new unsigned long[mSize];
       unsigned long * lFirstSiteIdForEachProc = new unsigned long[mSize];
       unsigned long * lNumberSitesPerProc = new unsigned long[mSize];


       unsigned long lTotalNonSolidSites = 0;

       for (unsigned int ii = 0; ii < bGlobalLatticeData.GetBlockCount(); ii++)
       {
       lTotalNonSolidSites += lNonSolidSitesPerBlock[ii];
       }

       GetInitialSiteDistribution(lFirstBlockIdForEachProc, lFirstSiteIdForEachProc,
       lNumberSitesPerProc, lTotalNonSolidSites, lNonSolidSitesPerBlock,
       bGlobalLatticeData);*/

      //TODO   ReadInBlocks

      MPI_File_close(&lFile);

      //OptimiseDomainDecomposition();

      double lMiddle = MPI_Wtime();

      DecomposeDomain(totalFluidSites, iReserveSteeringCore, bNetTop, bGlobalLatticeData);

      double lEnd = MPI_Wtime();

      *oReadTime = lMiddle - lStart;
      *oDecomposeTime = lEnd - lMiddle;
    }

    /**
     * Read in ALL of the blocks in the file.
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
    void TopologyReader::ReadAllBlocks(lb::GlobalLatticeData* bGlobLatDat,
                                       const unsigned int* bytesPerBlock,
                                       int &totalFluidSites,
                                       unsigned int siteMins[3],
                                       unsigned int siteMaxes[3],
                                       MPI_File iFile)
    {
      // Initialise counter variables.
      totalFluidSites = 0;

      siteMins[0] = UINT_MAX;
      siteMins[1] = UINT_MAX;
      siteMins[2] = UINT_MAX;
      siteMaxes[0] = 0;
      siteMaxes[1] = 0;
      siteMaxes[2] = 0;

      // Each site can have at most
      // * an unsigned int (config)
      // * 4 doubles
      // * 14 further doubles (in D3Q15)
      unsigned int maxBytesPerBlock = (bGlobLatDat->SitesPerBlockVolumeUnit) * (4 * 1 + 8 * (4
          + D3Q15::NUMVECTORS - 1));
      const unsigned int BlocksToReadInOneGo = 10;
      char* readBuffer = new char[maxBytesPerBlock * BlocksToReadInOneGo];

      // This makes sure we do the right number of read operations.
      for (unsigned int readNum = 0; readNum
          <= (bGlobLatDat->GetBlockCount() / BlocksToReadInOneGo); ++readNum)
      {
        const unsigned int upperLimitBlockNumber =
            util::NumericalFunctions::min<unsigned int>(bGlobLatDat->GetBlockCount(), (readNum + 1)
                * BlocksToReadInOneGo);
        const unsigned int lowerLimitBlockNumber = readNum * BlocksToReadInOneGo;

        unsigned int bytesToRead = 0;
        for (unsigned int ii = lowerLimitBlockNumber; ii < upperLimitBlockNumber; ++ii)
        {
          bytesToRead += bytesPerBlock[ii];
        }

        MPI_Status lStatus;

        MPI_File_read_all(iFile, readBuffer, (int) bytesToRead, MPI_BYTE, &lStatus);

        io::XdrMemReader lReader(readBuffer, bytesToRead);

        for (lb::BlockCounter lBlock(bGlobLatDat, lowerLimitBlockNumber); lBlock
            < upperLimitBlockNumber; lBlock++)
        {
          bGlobLatDat->Blocks[lBlock].site_data = NULL;
          bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite = NULL;
          bGlobLatDat->Blocks[lBlock].wall_data = NULL;

          if (bytesPerBlock[lBlock] == 0)
            continue;
          // Block contains some non-solid sites

          bGlobLatDat->Blocks[lBlock].site_data
              = new unsigned int[bGlobLatDat->SitesPerBlockVolumeUnit];
          bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite
              = new int[bGlobLatDat->SitesPerBlockVolumeUnit];

          int m = -1;

          for (unsigned int ii = 0; ii < bGlobLatDat->GetBlockSize(); ii++)
          {
            unsigned int site_i = lBlock.GetICoord(ii);

            for (unsigned int jj = 0; jj < bGlobLatDat->GetBlockSize(); jj++)
            {
              unsigned int site_j = lBlock.GetJCoord(jj);

              for (unsigned int kk = 0; kk < bGlobLatDat->GetBlockSize(); kk++)
              {
                unsigned int site_k = lBlock.GetKCoord(kk);

                ++m;

                unsigned int *site_type = &bGlobLatDat->Blocks[lBlock].site_data[m];
                lReader.readUnsignedInt(*site_type);

                if ( (*site_type & SITE_TYPE_MASK) == hemelb::lb::SOLID_TYPE)
                {
                  bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[m] = 1 << 30;
                  continue;
                }
                bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[m] = -1;

                ++totalFluidSites;

                siteMins[0] = hemelb::util::NumericalFunctions::min<unsigned int>(siteMins[0],
                                                                                  site_i);
                siteMins[1] = hemelb::util::NumericalFunctions::min<unsigned int>(siteMins[1],
                                                                                  site_j);
                siteMins[2] = hemelb::util::NumericalFunctions::min<unsigned int>(siteMins[2],
                                                                                  site_k);
                siteMaxes[0] = hemelb::util::NumericalFunctions::max<unsigned int>(siteMaxes[0],
                                                                                   site_i);
                siteMaxes[1] = hemelb::util::NumericalFunctions::max<unsigned int>(siteMaxes[1],
                                                                                   site_j);
                siteMaxes[2] = hemelb::util::NumericalFunctions::max<unsigned int>(siteMaxes[2],
                                                                                   site_k);

                if (bGlobLatDat->GetCollisionType(*site_type) != FLUID)
                {
                  // Neither solid nor simple fluid
                  if (bGlobLatDat->Blocks[lBlock].wall_data == NULL)
                  {
                    bGlobLatDat->Blocks[lBlock].wall_data
                        = new hemelb::lb::WallData[bGlobLatDat->SitesPerBlockVolumeUnit];
                  }

                  if (bGlobLatDat->GetCollisionType(*site_type) & INLET
                      || bGlobLatDat->GetCollisionType(*site_type) & OUTLET)
                  {
                    double temp;
                    // INLET or OUTLET or both.
                    // These values are the boundary normal and the boundary distance.
                    for (int l = 0; l < 3; l++)
                      lReader.readDouble(temp);

                    lReader.readDouble(temp);
                  }

                  if (bGlobLatDat->GetCollisionType(*site_type) & EDGE)
                  {
                    // EDGE bit set
                    for (int l = 0; l < 3; l++)
                      lReader.readDouble(bGlobLatDat->Blocks[lBlock].wall_data[m].wall_nor[l]);

                    double temp;
                    lReader.readDouble(temp);
                  }

                  for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
                    lReader.readDouble(bGlobLatDat->Blocks[lBlock].wall_data[m].cut_dist[l]);
                }
              } // kk
            } // jj
          } // ii

        }
      }

      delete[] readBuffer;
    }

    void TopologyReader::DecomposeDomain(int iTotalFluidSites,
                                         bool iReserveSteeringCore,
                                         NetworkTopology* bNetTop,
                                         const lb::GlobalLatticeData* bGlobLatDat)
    {
      // Allocations.  fluid sites will store actual number of fluid
      // sites per proc.  Site location will store up to 10000 of some
      // sort of coordinate.

      // Initialise the count of fluid sites on each processor to 0.
      bNetTop->FluidSitesOnEachProcessor = new int[bNetTop->GetProcessorCount()];

      for (unsigned int n = 0; n < bNetTop->GetProcessorCount(); n++)
      {
        bNetTop->FluidSitesOnEachProcessor[n] = 0;
      }

      // Count of fluid sites not yet visited.
      int lUnvisitedFluidSiteCount = iTotalFluidSites;

      // If one machine or one machine per proc.
      if (bNetTop->GetMachineCount() == 1 || bNetTop->GetMachineCount()
          == bNetTop->GetProcessorCount())
      {
        // Fluid sites per rank.
        int fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
            / (double) bNetTop->GetProcessorCount());

        //Rank we're looking at.
        int proc_count = 0;

        // If we're steering with more than one processor, save one processor for doing that.
        if (iReserveSteeringCore && bNetTop->GetProcessorCount() != 1)
        {
          fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
              / (double) (bNetTop->GetProcessorCount() - 1));
          proc_count = 1;
        }

        // In the simple case, simply divide fluid sites up between processors.
        AssignFluidSitesToProcessors(proc_count, fluid_sites_per_unit, lUnvisitedFluidSiteCount,
                                     -1, false, bNetTop, bGlobLatDat);
      }
      else
      {
        // Rank we are looking at.
        int proc_count = bNetTop->GetProcessorCount();

        // The number of fluid sites per machine is found by simple division, multiplied
        // by a factor to compensate for rank 0 being used for other things.
        double weight = (double) (bNetTop->GetProcessorCount())
            / (double) (bNetTop->GetProcessorCount() - 1);
        // Fluid sites per rank.
        int fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites * weight
            / bNetTop->GetMachineCount());

        // First, divide the sites up between machines.
        AssignFluidSitesToProcessors(proc_count, fluid_sites_per_unit, lUnvisitedFluidSiteCount,
                                     -1, true, bNetTop, bGlobLatDat);

        fluid_sites_per_unit = (int) ceil((double) lUnvisitedFluidSiteCount
            / (double) (bNetTop->GetProcessorCount() - 1));
        proc_count = 1;

        // For each machine, divide up the sites it has between its cores.
        for (unsigned int lMachineNumber = 0; lMachineNumber < bNetTop->GetMachineCount(); lMachineNumber++)
        {
          AssignFluidSitesToProcessors(proc_count, fluid_sites_per_unit, lUnvisitedFluidSiteCount,
                                       bNetTop->GetProcessorCount() + lMachineNumber, false,
                                       bNetTop, bGlobLatDat);
        }
      }

    }

    void TopologyReader::AssignFluidSitesToProcessors(int & proc_count,
                                                      int & iSitesPerProc,
                                                      int & bUnassignedSites,
                                                      const int iMarker,
                                                      const bool iIsMachineLevel,
                                                      NetworkTopology* bNetTop,
                                                      const lb::GlobalLatticeData* iGlobLatDat)
    {
      std::vector<SiteLocation*> *lSiteLocationA = new std::vector<SiteLocation*>;
      std::vector<SiteLocation*> *lSiteLocationB = new std::vector<SiteLocation*>;

      int lBlockNumber = -1;

      // Domain Decomposition.  Pick a site. Set it to the rank we are
      // looking at. Find its neighbours and put those on the same
      // rank, then find the next-nearest neighbours, etc. until we
      // have a completely joined region, or there are enough fluid
      // sites on the rank.  In the former case, start again at
      // another site. In the latter case, move on to the next rank.
      // Do this until all sites are assigned to a rank. There is a
      // high chance of of all sites on a rank being joined.

      int lSitesOnCurrentProc = 0;

      // Iterate over all blocks.
      for (unsigned int lBlockCoordI = 0; lBlockCoordI < iGlobLatDat->GetXBlockCount(); lBlockCoordI++)
      {
        for (unsigned int lBlockCoordJ = 0; lBlockCoordJ < iGlobLatDat->GetYBlockCount(); lBlockCoordJ++)
        {
          for (unsigned int lBlockCoordK = 0; lBlockCoordK < iGlobLatDat->GetZBlockCount(); lBlockCoordK++)
          {
            // Block number is the number of the block we're currently on.
            lBlockNumber++;

            // Point to a block of ProcessorRankForEachBlockSite.  If we are in a block of solids, move on.
            int *lProcRankForSite = iGlobLatDat->Blocks[lBlockNumber].ProcessorRankForEachBlockSite;

            // If the array of proc rank for each site is NULL, we're on an all-solid block.
            if (lProcRankForSite == NULL)
            {
              continue;
            }

            // Create variables for the index of this site on the block and the number of fluid sites
            // that have been assigned to the current processor.
            int lSiteNumber = -1;

            // For each dimension of the site co-ordinates, iterate over all values of the site
            // co-ordinates on the current block.
            for (unsigned int lSiteCoordI = lBlockCoordI * iGlobLatDat->GetBlockSize(); lSiteCoordI
                < lBlockCoordI * iGlobLatDat->GetBlockSize() + iGlobLatDat->GetBlockSize(); lSiteCoordI++)
            {
              for (unsigned int lSiteCoordJ = lBlockCoordJ * iGlobLatDat->GetBlockSize(); lSiteCoordJ
                  < lBlockCoordJ * iGlobLatDat->GetBlockSize() + iGlobLatDat->GetBlockSize(); lSiteCoordJ++)
              {
                for (unsigned int lSiteCoordK = lBlockCoordK * iGlobLatDat->GetBlockSize(); lSiteCoordK
                    < lBlockCoordK * iGlobLatDat->GetBlockSize() + iGlobLatDat->GetBlockSize(); lSiteCoordK++)
                {
                  // Keep track of the site number.
                  lSiteNumber++;

                  // Move on if the site is solid (ProcessorRankForEachBlockSite = BIG_NUMBER2) or has
                  // already been assigned to a rank (0 <= ProcessorRankForEachBlockSite < BIG_NUMBER2).
                  if (lProcRankForSite[lSiteNumber] != iMarker)
                  {
                    continue;
                  }
                  // We have found an unvisited fluid site to start growing the subdomain from.
                  // Assign it to the rank and update the fluid site counters.
                  lProcRankForSite[lSiteNumber] = proc_count;

                  ++lSitesOnCurrentProc;

                  // Record the location of this initial site.
                  lSiteLocationA->clear();
                  SiteLocation *lNew = new SiteLocation();
                  lNew->i = lSiteCoordI;
                  lNew->j = lSiteCoordJ;
                  lNew->k = lSiteCoordK;
                  lSiteLocationA->push_back(lNew);

                  // The subdomain can grow.
                  bool lIsRegionGrowing = true;

                  // While the region can grow (i.e. it is not bounded by solids or visited
                  // sites), and we need more sites on this particular rank.
                  while (lSitesOnCurrentProc < iSitesPerProc && lIsRegionGrowing)
                  {
                    for (unsigned int ii = 0; ii < lSiteLocationB->size(); ii++)
                    {
                      delete lSiteLocationB->operator [](ii);
                    }
                    lSiteLocationB->clear();

                    // Sites added to the edge of the mClusters during the iteration.
                    lIsRegionGrowing = false;

                    // For sites on the edge of the domain (sites_a), deal with the neighbours.
                    for (unsigned int index_a = 0; index_a < lSiteLocationA->size()
                        && lSitesOnCurrentProc < iSitesPerProc; index_a++)
                    {
                      lNew = lSiteLocationA->operator [](index_a);

                      for (unsigned int l = 1; l < D3Q15::NUMVECTORS && lSitesOnCurrentProc
                          < iSitesPerProc; l++)
                      {
                        // Record neighbour location.
                        int neigh_i = lNew->i + D3Q15::CX[l];
                        int neigh_j = lNew->j + D3Q15::CY[l];
                        int neigh_k = lNew->k + D3Q15::CZ[l];

                        // Move on if neighbour is outside the bounding box.
                        if (neigh_i == -1 || neigh_i == (int) iGlobLatDat->GetXSiteCount())
                          continue;
                        if (neigh_j == -1 || neigh_j == (int) iGlobLatDat->GetYSiteCount())
                          continue;
                        if (neigh_k == -1 || neigh_k == (int) iGlobLatDat->GetZSiteCount())
                          continue;

                        // Move on if the neighbour is in a block of solids (in which case
                        // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid or has already
                        // been assigned to a rank (in which case ProcessorRankForEachBlockSite != -1).  ProcessorRankForEachBlockSite
                        // was initialized in lbmReadConfig in io.cc.

                        // Pointer to the rank on which a particular fluid site
                        // resides.
                        int * proc_id_p = iGlobLatDat->GetProcIdFromGlobalCoords(neigh_i, neigh_j,
                                                                                 neigh_k);

                        if (proc_id_p == NULL || *proc_id_p != iMarker)
                        {
                          continue;
                        }
                        // Set the rank for a neighbour and update the fluid site counters.
                        *proc_id_p = proc_count;
                        lSitesOnCurrentProc++;

                        // Neighbour was found, so the region can grow.
                        lIsRegionGrowing = true;

                        // Record the location of the neighbour.
                        SiteLocation * lNewB = new SiteLocation();
                        lNewB->i = neigh_i;
                        lNewB->j = neigh_j;
                        lNewB->k = neigh_k;
                        lSiteLocationB->push_back(lNewB);
                      }
                    }
                    // When the new layer of edge sites has been found, swap the buffers for
                    // the current and new layers of edge sites.
                    std::vector<SiteLocation*> *tempP = lSiteLocationA;
                    lSiteLocationA = lSiteLocationB;
                    lSiteLocationB = tempP;
                  }

                  // If we have enough sites, we have finished.
                  if (lSitesOnCurrentProc >= iSitesPerProc)
                  {
                    if (!iIsMachineLevel)
                    {
                      bNetTop->FluidSitesOnEachProcessor[proc_count] = lSitesOnCurrentProc;
                    }

                    ++proc_count;
                    if (!iIsMachineLevel)
                    {
                      bUnassignedSites -= lSitesOnCurrentProc;
                      iSitesPerProc = (int) ceil((double) bUnassignedSites
                          / (double) (bNetTop->GetProcessorCount() - proc_count));
                    }
                    lSitesOnCurrentProc = 0;
                  }
                  // If not, we have to start growing a different region for the same rank:
                  // region expansions could get trapped.

                } // Site co-ord k
              } // Site co-ord j
            } // Site co-ord i
          } // Block co-ord k
        } // Block co-ord j
      } // Block co-ord i

      for (unsigned int ii = 0; ii < lSiteLocationA->size(); ii++)
      {
        delete lSiteLocationA->operator [](ii);
      }
      for (unsigned int ii = 0; ii < lSiteLocationB->size(); ii++)
      {
        delete lSiteLocationB->operator [](ii);
      }

      delete lSiteLocationA;
      delete lSiteLocationB;
    }

  }
}
