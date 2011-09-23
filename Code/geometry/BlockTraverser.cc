#include "geometry/BlockTraverser.h"
#include "geometry/LatticeData.h"
#include "util/Vector3D.h"
#include "vis/rayTracer/RayTracer.h"

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

    site_t BlockTraverser::CurrentBlockNumber()
    {
      return GetCurrentIndex();
    }

    util::Vector3D<site_t> BlockTraverser::GetSiteCoordinatesOfLowestSiteInCurrentBlock()
    {
      return GetCurrentLocation() * mLatticeData.GetBlockSize();
    }

    geometry::LatticeData::BlockData *
    BlockTraverser::GetCurrentBlockData()
    {
      return mLatticeData.GetBlock(mCurrentNumber);
    }

    geometry::LatticeData::BlockData *
    BlockTraverser::GetBlockDataForLocation(const util::Vector3D<site_t>& iLocation)
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
      return mLatticeData.IsValidBlockSite(iBlock.x, iBlock.y, iBlock.z);
    }

    bool BlockTraverser::GoToNextBlock()
    {
      return TraverseOne();
    }

    site_t BlockTraverser::GetXCount()
    {
      return mLatticeData.GetXBlockCount();
    }

    site_t BlockTraverser::GetYCount()
    {
      return mLatticeData.GetYBlockCount();
    }

    site_t BlockTraverser::GetZCount()
    {
      return mLatticeData.GetZBlockCount();
    }
  }
}
