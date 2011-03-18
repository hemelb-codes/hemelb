#ifndef HEMELB_TOPOLOGY_TOPOLOGYREADER_H
#define HEMELB_TOPOLOGY_TOPOLOGYREADER_H

#include "mpiInclude.h"
#include "lb/LbmParameters.h"
#include "geometry/GlobalLatticeData.h"
#include "net/net.h"

namespace hemelb
{

  namespace topology
  {

    class TopologyReader
    {
      public:
        TopologyReader(const bool reserveSteeringCore);
        ~TopologyReader();

        void LoadAndDecompose(geometry::GlobalLatticeData* bGlobalLatticeData,
                              int &totalFluidSites,
                              unsigned int siteMins[3],
                              unsigned int siteMaxes[3],
                              NetworkTopology* bNetTop,
                              lb::LbmParameters* bLbmParams,
                              SimConfig* bSimConfig,
                              double* lReadTime,
                              double* lDecomposeTime);

      private:
        struct BlockLocation
        {
            short int i, j, k;
        };

        void ReadPreamble(MPI_File xiFile,
                          lb::LbmParameters* bParams,
                          geometry::GlobalLatticeData* bGlobalLatticeData);

        void ReadHeader(MPI_File xiFile,
                        unsigned int iBlockCount,
                        unsigned int* sitesInEachBlock,
                        unsigned int* bytesUsedByBlockInDataFile);

        void BlockDecomposition(const unsigned int iBlockCount,
                                const geometry::GlobalLatticeData* iGlobLatDat,
                                const unsigned int* fluidSitePerBlock,
                                int* initialProcForEachBlock);

        void DivideBlocks(unsigned int unassignedBlocks,
                          unsigned int totalBlockCount,
                          unsigned int unitCount,
                          unsigned int* blocksOnEachUnit,
                          int* unitForEachBlock,
                          const unsigned int* fluidSitesPerBlock,
                          const geometry::GlobalLatticeData* iGlobLatDat);

        void ReadInLocalBlocks(MPI_File iFile,
                               const unsigned int* bytesPerBlock,
                               const int* unitForEachBlock,
                               const unsigned int localRank,
                               const geometry::GlobalLatticeData* iGlobLatDat);

        void OptimiseDomainDecomposition(const unsigned int* sitesPerBlock,
                                         const int* procForEachBlock,
                                         SimConfig* bSimConfig,
                                         lb::LbmParameters* bLbmParams,
                                         geometry::GlobalLatticeData* bGlobLatDat);

        MPI_Comm mTopologyComm;
        MPI_Group mTopologyGroup;
        int mTopologyRank;
        unsigned int mTopologySize;
        int mGlobalRank;
    };

  }

}

#endif /* HEMELB_TOPOLOGY_TOPOLOGYREADER_H */
