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
    */
    class Decomposition
    {
    public:
      Decomposition(const site_t BlockCount, bool *readBlock, const proc_t readingGroupSize, 
        const topology::NetworkTopology &atopology=*topology::NetworkTopology::Instance() ); // Temporarily during the refactor, constructed just to abstract the block sharing bit
      std::vector<proc_t> ProcessorsNeedingBlock(const site_t &block){return procsWantingBlocksBuffer[block];}
    private:
      std::vector<std::vector<proc_t> > procsWantingBlocksBuffer;
      net::Net net;
      const topology::NetworkTopology &topology;
    };
  }
}
#endif // HEMELB_GEOMETRY_DECOMPOSITION_H