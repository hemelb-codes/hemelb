#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis 
  {
    namespace raytracer 
    {
      SiteTraverser::SiteTraverser(const geometry::LatticeData::LatticeData * iLatticeData,
							    const site_t iBlockId)
	: mLatticeData(iLatticeData),
	  mBlockId(iBlockId)
      { 	
      }
      
      site_t SiteTraverser::GetXCount()
      {
	return GetBlockSize();
      }

      site_t SiteTraverser::GetYCount()
      {
	return GetBlockSize();
      }

      site_t SiteTraverser::GetZCount()
      {
	return GetBlockSize();
      }

      site_t SiteTraverser::GetBlockSize()
      {
	return mLatticeData->GetBlockSize();
      }	
    }
  }
}
