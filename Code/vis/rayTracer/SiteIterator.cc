#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis 
  {
    namespace raytracer 
    {
      RayTracer::ClusterBuilder::SiteIterator::SiteIterator(const geometry::LatticeData::LatticeData * iLatticeData,
							    const site_t iBlockId)
	: mLatticeData(iLatticeData),
	  mBlockId(iBlockId)
      { 	
      }
      

      site_t RayTracer::ClusterBuilder::SiteIterator::GetXCount()
      {
	return GetBlockSize();
      }

      site_t RayTracer::ClusterBuilder::SiteIterator::GetYCount()
      {
	return GetBlockSize();
      }

      site_t RayTracer::ClusterBuilder::SiteIterator::GetZCount()
      {
	return GetBlockSize();
      }

      site_t RayTracer::ClusterBuilder::SiteIterator::GetBlockSize()
      {
	return mLatticeData->GetBlockSize();
      }	
    }
  }
}
