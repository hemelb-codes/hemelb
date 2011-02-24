#ifndef HEMELB_TOPOLOGY_TOPOLOGYREADER_H
#define HEMELB_TOPOLOGY_TOPOLOGYREADER_H

#include "mpiInclude.h"
#include "lb/LbmParameters.h"
#include "lb/GlobalLatticeData.h"
#include "net/net.h"

namespace hemelb
{

  namespace topology
  {

    class TopologyReader
    {
      public:
        TopologyReader();
        ~TopologyReader();

        void LoadAndDecompose(lb::GlobalLatticeData &bGlobalLatticeData,
                              net::Net *net,
                              lb::LbmParameters * bLbmParams,
                              SimConfig * bSimConfig);

        void ReadPreamble(MPI_File xiFile,
                          lb::LbmParameters * bParams,
                          lb::GlobalLatticeData &bGlobalLatticeData);

        void ReadHeader(MPI_File xiFile,
                        unsigned int iBlockCount,
                        unsigned int* sitesInEachBlock,
                        unsigned int* bytesUsedByBlockInDataFile);

      private:

        /*   void
         GetNeighbourLocation(int lFromBlockId,
         int lFromSiteId,
         int lDirection,
         hemelb::lb::GlobalLatticeData &bGlobalLatticeData,
         int * lToBlockId,
         int * lToSiteId);*/

        void GetNonSolidSitesPerBlock(int bNonSolidSitesPerBlock[],
                                      net::Net *iNet,
                                      MPI_File iFile,
                                      const lb::GlobalLatticeData &bGlobalLatticeData);
        void GetInitialSiteDistribution(unsigned long oFirstBlockIdForEachProc[],
                                        unsigned long oFirstSiteIdForEachProc[],
                                        unsigned long oNumberSitesPerProc[],
                                        unsigned long oTotalSiteCount,
                                        const int iNonSolidSitesPerBlock[],
                                        const hemelb::lb::GlobalLatticeData & iGlobLatDat);

        void ReadInBlocks(const unsigned long iFirstBlockIdForEachProc[],
                          const unsigned long iFirstSiteNumberForEachProc[],
                          const unsigned long iNumberSitesPerProc[],
                          const unsigned long iTotalSiteCount,
                          const int iNonSolidSitesPerBlock[],
                          const hemelb::lb::GlobalLatticeData & iGlobLatDat);

        void OptimiseDomainDecomposition(const unsigned long iFirstBlockIdForEachProc[],
                                         const unsigned long iFirstSiteIdForEachProc[],
                                         const unsigned long iNumberSitesPerProc[],
                                         const unsigned long iTotalSiteCount,
                                         const int iNonSolidSitesPerBlock[],
                                         const hemelb::lb::GlobalLatticeData & iGlobLatDat);

        MPI_Comm mCommunicator;
        MPI_Group mGroup;
        int mRank;
        int mSize;
    };

  }

}

#endif /* HEMELB_TOPOLOGY_TOPOLOGYREADER_H */
