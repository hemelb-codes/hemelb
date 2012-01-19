#ifndef HEMELB_GEOMETRY_BLOCKTRAVERSER_H
#define HEMELB_GEOMETRY_BLOCKTRAVERSER_H

#include "D3Q15.h"
#include "geometry/Block.h"
#include "geometry/VolumeTraverser.h"
#include "geometry/SiteTraverser.h"

namespace hemelb
{
  namespace geometry
  {
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

        const Block& GetCurrentBlockData();

        const Block& GetBlockDataForLocation(const util::Vector3D<site_t>& iLocation);

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
