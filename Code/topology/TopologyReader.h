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

        void LoadAndDecompose(lb::GlobalLatticeData* bGlobalLatticeData,
                              int &totalFluidSites,
                              unsigned int siteMins[3],
                              unsigned int siteMaxes[3],
                              lb::LbmParameters* bLbmParams,
                              SimConfig* bSimConfig);

        void DecomposeDomain(int iTotalFluidSites,
                             bool iReserveSteeringCore,
                             NetworkTopology* bNetTop,
                             const lb::GlobalLatticeData & bGlobLatDat);

      private:
        struct SiteLocation
        {
            short int i, j, k;
        };

        void AssignFluidSitesToProcessors(int & proc_count,
                                          int & fluid_sites_per_unit,
                                          int & unvisited_fluid_sites,
                                          const int iCurrentProcId,
                                          const bool iIsMachineLevel,
                                          NetworkTopology* bNetTop,
                                          const lb::GlobalLatticeData &iGlobLatDat);

        void ReadPreamble(MPI_File xiFile,
                          lb::LbmParameters* bParams,
                          lb::GlobalLatticeData* bGlobalLatticeData);

        void ReadHeader(MPI_File xiFile,
                        unsigned int iBlockCount,
                        unsigned int* sitesInEachBlock,
                        unsigned int* bytesUsedByBlockInDataFile);

        void ReadAllBlocks(lb::GlobalLatticeData* bGlobLatDat,
                           const unsigned int* bytesPerBlock,
                           int &totalFluidSites,
                           unsigned int siteMins[3],
                           unsigned int siteMaxes[3],
                           MPI_File iFile);

        /*   void
         GetNeighbourLocation(int lFromBlockId,
         int lFromSiteId,
         int lDirection,
         hemelb::lb::GlobalLatticeData &bGlobalLatticeData,
         int * lToBlockId,
         int * lToSiteId);*/

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
