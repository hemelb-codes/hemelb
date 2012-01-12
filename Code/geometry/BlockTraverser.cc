#include "geometry/BlockTraverser.h"
#include "geometry/LatticeData.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    BlockTraverser::BlockTraverser(const geometry::LatticeData& iLatDat) :
        VolumeTraverser(), mLatticeData(iLatDat)
    {
    }

    BlockTraverser::~BlockTraverser()
    {
    }

    site_t BlockTraverser::CurrentBlockNumber() const
    {
      return GetCurrentIndex();
    }

    util::Vector3D<site_t> BlockTraverser::GetSiteCoordinatesOfLowestSiteInCurrentBlock()
    {
      return GetCurrentLocation() * mLatticeData.GetBlockSize();
    }

    const geometry::BlockData* BlockTraverser::GetCurrentBlockData()
    {
      return mLatticeData.GetBlock(GetCurrentIndex());
    }

    const geometry::BlockData* BlockTraverser::GetBlockDataForLocation(const util::Vector3D<site_t>& iLocation)
    {
      return mLatticeData.GetBlock(GetIndexFromLocation(iLocation));
    }

    site_t BlockTraverser::GetBlockSize()
    {
      return mLatticeData.GetBlockSize();
    }

    SiteTraverser BlockTraverser::GetSiteTraverser()
    {
      return SiteTraverser(mLatticeData);
    }

    bool BlockTraverser::IsValidLocation(util::Vector3D<site_t> iBlock)
    {
      return mLatticeData.IsValidBlock(iBlock.x, iBlock.y, iBlock.z);
    }

    bool BlockTraverser::GoToNextBlock()
    {
      return TraverseOne();
    }

    site_t BlockTraverser::GetXCount() const
    {
      return mLatticeData.GetXBlockCount();
    }

    site_t BlockTraverser::GetYCount() const
    {
      return mLatticeData.GetYBlockCount();
    }

    site_t BlockTraverser::GetZCount() const
    {
      return mLatticeData.GetZBlockCount();
    }
  }
}
