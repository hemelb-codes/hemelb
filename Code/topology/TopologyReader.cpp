#include "TopologyReader.h"
#include "io/XdrMemReader.h"
extern "C"
{
#include "parmetis/parmetislib.h"
}

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
    void TopologyReader::ReadPreamble(MPI_File xiFile,
                                      hemelb::lb::LbmParameters * bParams,
                                      hemelb::lb::GlobalLatticeData &bGlobalLatticeData)
    {
      std::string lMode = "native";

      MPI_File_set_view(xiFile, 0, MPI_BYTE, MPI_BYTE, &lMode[0], MPI_INFO_NULL);


      // OLD Version

      static const int mPreambleBytes = 24;

      char lPreambleBuffer[mPreambleBytes];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lPreambleBuffer, mPreambleBytes, MPI_BYTE, &lStatus);
      hemelb::io::XdrReader myReader = hemelb::io::XdrMemReader(lPreambleBuffer, mPreambleBytes);
      // Not the ideal way to do this, but has to be this way as the old system used
      // doubles for the stress type. -1.0 signified shear stress, 1.0 meant von Mises.
      double lTempStressType;

      myReader.readDouble(lTempStressType);
      bParams->StressType = (lTempStressType == -1.0)
        ? hemelb::lb::ShearStress
        : ( (lTempStressType == 1.0)
          ? hemelb::lb::VonMises
          : hemelb::lb::IgnoreStress);

      int lBlocksX, lBlocksY, lBlocksZ, lBlockSize;

      myReader.readInt(lBlocksX);
      myReader.readInt(lBlocksY);
      myReader.readInt(lBlocksZ);
      myReader.readInt(lBlockSize);

      bGlobalLatticeData.SetBasicDetails(lBlocksX, lBlocksY, lBlocksZ, lBlockSize);

      // New version.
      //      // The config file starts with:
      //      // * 1 unsigned int for stress type
      //      // * 3 unsigned ints for the number of blocks in the x, y, z directions
      //      // * 1 unsigned int for the block size (number of sites along one edge of a block)
      //      // * 1 unsigned int for the total number of blocks in the cube, which should be equal to the
      //      //   product of the x, y and z block counts.
      //      static const int PreambleBytes = 24;
      //
      //      char lPreambleBuffer[PreambleBytes];
      //
      //      MPI_Status lStatus;
      //
      //      MPI_File_read_all(xiFile, lPreambleBuffer, PreambleBytes, MPI_BYTE, &lStatus);
      //
      //      hemelb::io::XdrReader preambleReader = hemelb::io::XdrMemReader(lPreambleBuffer,
      //                                                                      PreambleBytes);
      //
      //      // Variables we'll read.
      //      unsigned int stressType, blocksX, blocksY, blocksZ, blockSize;
      //
      //      preambleReader.readUnsignedInt(stressType);
      //      preambleReader.readUnsignedInt(blocksX);
      //      preambleReader.readUnsignedInt(blocksY);
      //      preambleReader.readUnsignedInt(blocksZ);
      //      preambleReader.readUnsignedInt(blockSize);
      //
      //      bParams->StressType = stressType;
      //
      //      bGlobalLatticeData.SetBasicDetails(blocksX, blocksY, blocksZ, blockSize);
      //
      //      unsigned int blockCount;
      //
      //      preambleReader.readUnsignedInt(blockCount);
      //
      //      if (bGlobalLatticeData.GetBlockCount() != blockCount)
      //      {
      //        printf(
      //               "Input file may be corrupted: Recorded block count was %i but should have been %i.\n",
      //               blockCount, bGlobalLatticeData.GetBlockCount());
      //        exit(1);
      //      }
    }

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
    }

    void TopologyReader::GetNonSolidSitesPerBlock(int bNonSolidSitesPerBlock[],
                                                  net::Net* iNet,
                                                  MPI_File iFile,
                                                  const hemelb::lb::GlobalLatticeData &bGlobalLatticeData)
    {
      // Each block has an int flag, each site has at most an unsigned int, 8 doubles, and (Num-vectors - 1) doubles.
      unsigned int lLengthPerBlock = 4 + bGlobalLatticeData.SitesPerBlockVolumeUnit * (4 + 8 * 8
          + 8 * (D3Q15::NUMVECTORS - 1));

      // 500000 bytes is a reasonable amount to read at once (a significant read, but small enough
      // that it's barely noticeable in terms of memory usage.
      unsigned int lBytesToReadAtOnce = (lLengthPerBlock > 500000)
        ? lLengthPerBlock
        : 500000;

      char * lBlockDataBuffer = new char[lBytesToReadAtOnce];

      MPI_Status lStatus;

      MPI_File_read_all(iFile, lBlockDataBuffer, lBytesToReadAtOnce, MPI_BYTE, &lStatus);

      hemelb::io::XdrMemReader myReader = hemelb::io::XdrMemReader(lBlockDataBuffer,
                                                                   lBytesToReadAtOnce);

      for (unsigned int n = 0; n < bGlobalLatticeData.GetBlockCount(); n++)
      {
        // If we've read so far through the buffer that we might go over at the end of
        // this block, read in some more data.
        if ( (myReader.GetPosition() + lLengthPerBlock) < lBytesToReadAtOnce)
        {
          for (unsigned int ii = myReader.GetPosition(); ii < lBytesToReadAtOnce; ii++)
          {
            lBlockDataBuffer[ii - myReader.GetPosition()] = lBlockDataBuffer[ii];
          }

          MPI_File_read_all(iFile,
                            &lBlockDataBuffer[ (lBytesToReadAtOnce - myReader.GetPosition())],
                            myReader.GetPosition(), MPI_BYTE, &lStatus);

          myReader = hemelb::io::XdrMemReader(lBlockDataBuffer, lBytesToReadAtOnce);
        }

        int flag;

        myReader.readInt(flag);

        // If block contains all solid sites, move on.
        if (flag == 0)
        {
          continue;
        }

        for (unsigned int m = 0; m < bGlobalLatticeData.SitesPerBlockVolumeUnit; m++)
        {
          unsigned int site_type;
          myReader.readUnsignedInt(site_type);

          if ( (site_type & SITE_TYPE_MASK) == hemelb::lb::SOLID_TYPE)
          {
            continue;
          }

          ++bNonSolidSitesPerBlock[n];

          // Have to read in all the rest of the data for this site, to pass it.
          if (bGlobalLatticeData.GetCollisionType(site_type) != FLUID)
          {
            double temp;
            // Neither solid nor simple fluid
            if (bGlobalLatticeData.GetCollisionType(site_type) & INLET
                || bGlobalLatticeData.GetCollisionType(site_type) & OUTLET)
            {
              // INLET or OUTLET or both.
              // These values are the boundary normal and the boundary distance.
              for (int l = 0; l < 3; l++)
                myReader.readDouble(temp);

              myReader.readDouble(temp);
            }

            if (bGlobalLatticeData.GetCollisionType(site_type) & EDGE)
            {
              // EDGE bit set
              for (int l = 0; l < 3; l++)
                myReader.readDouble(temp);

              double temp;
              myReader.readDouble(temp);
            }

            for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
              myReader.readDouble(temp);
          }
        } // m
      } // n

      delete[] lBlockDataBuffer;
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

      //(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, floattype *tpwgts, floattype *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
      //ParMETIS_V3_PartKway();
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

              if (lPutativeBlock >= 0 && lPutativeBlock < (int)iGlobLatDat.GetBlockCount())
              {
                lReadBlock[lPutativeBlock] = true;
              }
            }
          }
        }
      }

      // Now we actually do the reading.

    }

    void TopologyReader::LoadAndDecompose(lb::GlobalLatticeData &bGlobalLatticeData,
                                          net::Net *net,
                                          lb::LbmParameters * bLbmParams,
                                          SimConfig * bSimConfig)
    {
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

      int * lNonSolidSitesPerBlock = new int[bGlobalLatticeData.GetBlockCount()];

      for (unsigned int ii = 0; ii < bGlobalLatticeData.GetBlockCount(); ii++)
      {
        // Set the non-solid sites to 0 for the current block.
        lNonSolidSitesPerBlock[ii] = 0;
      }

      GetNonSolidSitesPerBlock(lNonSolidSitesPerBlock, net, lFile, bGlobalLatticeData);

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
                                 bGlobalLatticeData);

      //TODO   ReadInBlocks


      //OptimiseDomainDecomposition();
    }

  }
}
