#ifndef HEMELB_TOPOLOGY_TOPOLOGYREADER_H
#define HEMELB_TOPOLOGY_TOPOLOGYREADER_H

#include "mpiInclude.h"
#include "lb/LbmParameters.h"
#include "lb/GlobalLatticeData.h"
#include "net.h"

namespace hemelb
{

  namespace topology
  {

    class TopologyReader
    {
      public:
        TopologyReader();
        ~TopologyReader();

        void
        LoadAndDecompose(hemelb::lb::GlobalLatticeData &bGlobalLatticeData,
                         Net *net,
                         hemelb::lb::LbmParameters * bLbmParams,
                         SimConfig * bSimConfig);

        void
        PreReadConfigFile(MPI_File xiFile,
                          hemelb::lb::LbmParameters * bParams,
                          hemelb::lb::GlobalLatticeData &bGlobalLatticeData);

      private:

        void
            GetNonSolidSitesPerBlock(int bNonSolidSitesPerBlock[],
                                     Net *iNet,
                                     MPI_File iFile,
                                     const hemelb::lb::GlobalLatticeData &bGlobalLatticeData);
        void
            GetInitialSiteDistribution(unsigned long oFirstBlockIdForEachProc[],
                                       unsigned long oFirstSiteIdForEachProc[],
                                       unsigned long oNumberSitesPerProc[],
                                       unsigned long oTotalSiteCount,
                                       const int iNonSolidSitesPerBlock[],
                                       const hemelb::lb::GlobalLatticeData & iGlobLatDat);

        void OptimiseDomainDecomposition();

        MPI_Comm mCommunicator;
        MPI_Group mGroup;
        int mRank;
        int mSize;

        // The config file starts with a double and 4 ints.
        // In Xdr this occupies 8 + 4 * 4 bytes.
        static const int mPreambleBytes = 24;

    };

  }

}

#endif /* HEMELB_TOPOLOGY_TOPOLOGYREADER_H */
