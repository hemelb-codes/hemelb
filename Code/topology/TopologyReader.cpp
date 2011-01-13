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

    void TopologyReader::PreReadConfigFile(MPI_File xiFile,
                                           hemelb::lb::LbmParameters * bParams,
                                           hemelb::lb::GlobalLatticeData &bGlobalLatticeData)
    {
      std::string lMode = "native";

      MPI_File_set_view(xiFile, 0, MPI_BYTE, MPI_BYTE, &lMode[0], MPI_INFO_NULL);

      char lPreambleBuffer[mPreambleBytes];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lPreambleBuffer, mPreambleBytes, MPI_BYTE,
                        &lStatus);

      hemelb::io::XdrReader myReader =
          hemelb::io::XdrMemReader(lPreambleBuffer, mPreambleBytes);

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

      bGlobalLatticeData.SetBasicDetails(lBlocksX, lBlocksY, lBlocksZ,
                                         lBlockSize);
    }

    void TopologyReader::GetNonSolidSitesPerBlock(int bNonSolidSitesPerBlock[],
                                                  Net* iNet,
                                                  MPI_File iFile,
                                                  const hemelb::lb::GlobalLatticeData &bGlobalLatticeData)
    {
      // Each block has an int flag, each site has at most an unsigned int, 8 doubles, and (Num-vectors - 1) doubles.
      unsigned int lLengthPerBlock = 4
          + bGlobalLatticeData.SitesPerBlockVolumeUnit * (4 + 8 * 8 + 8
              * (D3Q15::NUMVECTORS - 1));

      // 500000 bytes is a reasonable amount to read at once (a significant read, but small enough
      // that it's barely noticeable in terms of memory usage.
      unsigned int lBytesToReadAtOnce = (lLengthPerBlock > 500000)
        ? lLengthPerBlock
        : 500000;

      char * lBlockDataBuffer = new char[lBytesToReadAtOnce];

      MPI_Status lStatus;

      MPI_File_read_all(iFile, lBlockDataBuffer, lBytesToReadAtOnce, MPI_BYTE,
                        &lStatus);

      hemelb::io::XdrMemReader myReader =
          hemelb::io::XdrMemReader(lBlockDataBuffer, lBytesToReadAtOnce);

      for (int n = 0; n < bGlobalLatticeData.GetBlockCount(); n++)
      {
        // If we've read so far through the buffer that we might go over at the end of
        // this block, read in some more data.
        if ( (myReader.GetPosition() + lLengthPerBlock) < lBytesToReadAtOnce)
        {
          for (unsigned int ii = myReader.GetPosition(); ii
              < lBytesToReadAtOnce; ii++)
          {
            lBlockDataBuffer[ii - myReader.GetPosition()]
                = lBlockDataBuffer[ii];
          }

          MPI_File_read_all(iFile, &lBlockDataBuffer[ (lBytesToReadAtOnce
              - myReader.GetPosition())], myReader.GetPosition(), MPI_BYTE,
                            &lStatus);

          myReader = hemelb::io::XdrMemReader(lBlockDataBuffer,
                                              lBytesToReadAtOnce);
        }

        int flag;

        myReader.readInt(flag);

        // If block contains all solid sites, move on.
        if (flag == 0)
        {
          continue;
        }

        for (int m = 0; m < bGlobalLatticeData.SitesPerBlockVolumeUnit; m++)
        {
          unsigned int site_type;
          myReader.readUnsignedInt(site_type);

          if ( (site_type & SITE_TYPE_MASK) == hemelb::lb::SOLID_TYPE)
          {
            continue;
          }

          ++bNonSolidSitesPerBlock[n];

          // Have to read in all the rest of the data for this site, to pass it.
          if (iNet->GetCollisionType(site_type) != FLUID)
          {
            double temp;
            // Neither solid nor simple fluid
            if (iNet->GetCollisionType(site_type) & INLET
                || iNet->GetCollisionType(site_type) & OUTLET)
            {
              // INLET or OUTLET or both.
              // These values are the boundary normal and the boundary distance.
              for (int l = 0; l < 3; l++)
                myReader.readDouble(temp);

              myReader.readDouble(temp);
            }

            if (iNet->GetCollisionType(site_type) & EDGE)
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

      for (int ii = 0; ii < iGlobLatDat.GetBlockCount(); ii++)
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

        for (; (lSitesToPass > 0) && (lCurrentBlockId
            <= iGlobLatDat.GetBlockCount());)
        {
          int lSitesLeftOnCurrentBlock =
              iNonSolidSitesPerBlock[lCurrentBlockId] - lCurrentSiteId;

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

    void TopologyReader::OptimiseDomainDecomposition()
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

      //ParMETIS_V3_PartKway    ()
    }

    void TopologyReader::LoadAndDecompose(hemelb::lb::GlobalLatticeData &bGlobalLatticeData,
                                          Net *net,
                                          hemelb::lb::LbmParameters * bLbmParams,
                                          SimConfig * bSimConfig)
    {
      MPI_File lFile;
      int lError;

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      lError = MPI_File_open(MPI_COMM_WORLD, &bSimConfig->DataFilePath[0],
                             MPI_MODE_RDONLY, MPI_INFO_NULL, &lFile);

      if (lError != 0)
      {
        fprintf(stderr, "Unable to open file %s [rank %i], exiting\n",
                bSimConfig->DataFilePath.c_str(), mRank);
        fflush(0x0);
        exit(0x0);
      }
      else
      {
        fprintf(stderr, "Opened config file %s [rank %i]\n",
                bSimConfig->DataFilePath.c_str(), mRank);
      }
      fflush(NULL);

      PreReadConfigFile(lFile, bLbmParams, bGlobalLatticeData);

      int * lNonSolidSitesPerBlock =
          new int[bGlobalLatticeData.GetBlockCount()];

      for (int ii = 0; ii < bGlobalLatticeData.GetBlockCount(); ii++)
      {
        // Set the non-solid sites to 0 for the current block.
        lNonSolidSitesPerBlock[ii] = 0;
      }

      GetNonSolidSitesPerBlock(lNonSolidSitesPerBlock, net, lFile,
                               bGlobalLatticeData);

      unsigned long * lFirstBlockIdForEachProc = new unsigned long[mSize];
      unsigned long * lFirstSiteIdForEachProc = new unsigned long[mSize];
      unsigned long * lNumberSitesPerProc = new unsigned long[mSize];

      unsigned long lTotalNonSolidSites = 0;

      for (int ii = 0; ii < bGlobalLatticeData.GetBlockCount(); ii++)
      {
        lTotalNonSolidSites += lNonSolidSitesPerBlock[ii];
      }

      GetInitialSiteDistribution(lFirstBlockIdForEachProc,
                                 lFirstSiteIdForEachProc, lNumberSitesPerProc,
                                 lTotalNonSolidSites, lNonSolidSitesPerBlock,
                                 bGlobalLatticeData);

      OptimiseDomainDecomposition();
    }

  }
}
