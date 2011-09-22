#include "geometry/SiteTraverser.h"

namespace hemelb
{
  namespace geometry 
  {
    SiteTraverser::SiteTraverser(const geometry::LatticeData::LatticeData& iLatticeData)
      : mLatticeData(iLatticeData)
	  
    { 	
    }

    SiteTraverser::~SiteTraverser()
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
