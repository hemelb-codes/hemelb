#ifndef HEMELB_GEOMETRY_BLOCKTRAVERSER_H
#define HEMELB_GEOMETRY_BLOCKTRAVERSER_H

#include "D3Q15.h"
#include "geometry/VolumeTraverser.h"
#include "geometry/SiteTraverser.h"

namespace hemelb
{
  namespace geometry
  {

    // Data about each global block in the lattice,
    // site_data[] is an array containing individual lattice site data
    // within a global block.
    struct BlockData
    {
        BlockData()
        {
        }

        BlockData(site_t sitesPerBlock)
        {
          processorRankForEachBlockSite.resize(sitesPerBlock, BIG_NUMBER3);
          localContiguousIndex.resize(sitesPerBlock, BIG_NUMBER3);
        }

        ~BlockData()
        {
        }

        // An array of the ranks on which each lattice site within the block resides.
        std::vector<proc_t> processorRankForEachBlockSite;

        // The local index for each site on the block in the LocalLatticeData.
        std::vector<site_t> localContiguousIndex;
    };

    /**
     *BlockTraverser is used to traverse the blocks in a lattice sequentially.
     */
    class BlockTraverser : public VolumeTraverser
    {
      public:
        BlockTraverser(const LatticeData& iLatDat);
        ~BlockTraverser();

        site_t CurrentBlockNumber() const;

        util::Vector3D<site_t> GetSiteCoordinatesOfLowestSiteInCurrentBlock();

        const BlockData* GetCurrentBlockData();

        const BlockData* GetBlockDataForLocation(const util::Vector3D<site_t>& iLocation);

        site_t GetBlockSize();

        SiteTraverser GetSiteTraverser();

        site_t GetXCount() const;

        site_t GetYCount() const;

        site_t GetZCount() const;

        bool IsValidLocation(util::Vector3D<site_t> block);

      protected:
        bool GoToNextBlock();

        const LatticeData & mLatticeData;
    };
  }
}

#endif // HEMELB_GEOMETRY_BLOCKTRAVERSER_H
