#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis 
  {
    namespace raytracer 
    {
      RayTracer::ClusterBuilder::SiteIterator::SiteIterator(const geometry::LatticeData:LatticeData * iLatticeData, site_t iBlockd)
      {
	mLatticeData = iLatticeData;
	mBlockId = iBlockId;
      }
      

      virtual site_t RayTracer::ClusterBuilder::SiteIterator::GetXCount()
      {
	mLatticeData->getBlockSize();
      }

      virtual site_t RayTracer::ClusterBuilder::SiteIterator::GetYCount()
      {
	mLatticeData->getBlockSize();
      }
      virtual site_t RayTracer::ClusterBuilder::SiteIterator::GetZCount()
      {
	mLatticeData->getBlockSize();
      }
	
    }
  }
}
