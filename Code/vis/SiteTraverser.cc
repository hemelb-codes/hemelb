#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis 
  {
    SiteTraverser::SiteTraverser(const geometry::LatticeData::LatticeData& iLatticeData)
      : mLatticeData(iLatticeData)
	  
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
      return mLatticeData.GetBlockSize();
    }	
  }
}
