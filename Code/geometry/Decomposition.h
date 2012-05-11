#ifndef HEMELB_GEOMETRY_DECOMPOSITION_H
#define HEMELB_GEOMETRY_DECOMPOSITION_H
#include <vector>
#include "net/net.h"
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace geometry
  {
    /*
     JH Note: Class created to make the new decomposition communication strategy in #142 testable.
     Will not initially contain all decomposition-related code, will gradually add this.
     Eventually, will do the decomposition, and have clean interfaces to provide each site's owned and needed blocks during geometry reading.
     Site ownership during simulation will not be managed by this class, but by the lattice data as currently.
     */
    /***
     *  Class defining HemeLB domain decomposition
     Used by geometry reader to know where to send which blocks.
     @tparam Net Class implementing the Net Communication Protocol, used to share information.
     */
    template<class Net> class DecompositionBase
    {
      public:
        /***
         * Constructor for Decomposition, in the temporary state where not all decomposition-related code is here.
         * @param BlockCount Count of blocks
         * @param readBlock Temporary input - which cores need which blocks, as an array of booleans.
         * @param areadingGroupSize Number of cores to use for reading blocks
         * @param anet Instance of Net Communication Protocol -- typically just Net(comm)
         * @param comm Communicator used for decomposition topology
         * @param rank Rank in decomposition topology
         * @param size Size of decomposiiton topology
         */
        DecompositionBase(const site_t BlockCount,
                          const std::vector<bool>& readBlock,
                          const proc_t areadingGroupSize,
                          Net &anet,
                          MPI_Comm comm,
                          const proc_t rank,
                          const proc_t size); // Temporarily during the refactor, constructed just to abstract the block sharing bit

        /***
         * Which processors need a given block?
         * @param block Block number to query
         * @return Vector of ranks in the decomposition topology which need this block
         */
        const std::vector<proc_t> & ProcessorsNeedingBlock(const site_t &block) const
        {
          return procsWantingBlocksBuffer[block];
        }

        /***
         * Which core should be responsible for reading a given block? This core does not necessarily
         * require information about the block
         *
         * @param blockNumber Block number to query
         * @return Rank in the decomposition topology, for core which should read the block.
         */
        proc_t GetReadingCoreForBlock(const site_t blockNumber) const;
      private:
        std::vector<std::vector<proc_t> > procsWantingBlocksBuffer;
        Net &net;
        MPI_Comm decompositionCommunicator;
        const proc_t decompositionCommunicatorRank;
        const proc_t decompositionCommunicatorSize;
        const proc_t readingGroupSize;
    };
    typedef DecompositionBase<net::Net> Decomposition;
  }
}
#endif // HEMELB_GEOMETRY_DECOMPOSITION_H
