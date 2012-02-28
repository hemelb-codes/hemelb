#ifndef HEMELB_GEOMETRY_DECOMPOSITION_H
#define HEMELB_GEOMETRY_DECOMPOSITION_H
#include <vector>
#include "net/net.h"
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace geometry
  {
    /**
    Class defining HemeLB domain decomposition
    Used by geometry reader to know where to send which blocks.
    **/
    /*
    JH Note: Class created to make the new decomposition communication strategy in #142 testable.
    Will not initially contain all decomposition-related code, will gradually add this.
    Eventually, will do the decomposition, and have clean interfaces to provide each site's owned and needed blocks during geometry reading.
    Site ownership during simulation will not be managed by this class, but by the lattice data as currently.
    */
    template<class Net> class DecompositionBase
    {
    public:
      DecompositionBase(const site_t BlockCount, bool *readBlock, const proc_t readingGroupSize, Net &anet, MPI_Comm comm,
        const proc_t rank, const proc_t size); // Temporarily during the refactor, constructed just to abstract the block sharing bit
      std::vector<proc_t> ProcessorsNeedingBlock(const site_t &block){
        return procsWantingBlocksBuffer[block];
      }
    private:
      std::vector<std::vector<proc_t> > procsWantingBlocksBuffer;
      Net &net;
      MPI_Comm decompositionCommunicator;
      const proc_t decompositionCommunicatorRank;
      const proc_t decompositionCommunicatorSize;
    };
    typedef DecompositionBase<net::Net> Decomposition;
  }
}
#endif // HEMELB_GEOMETRY_DECOMPOSITION_H
