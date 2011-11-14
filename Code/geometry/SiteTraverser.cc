#include "geometry/LatticeData.h"
#include "geometry/SiteTraverser.h"

namespace hemelb
{
  namespace geometry
  {
    SiteTraverser::SiteTraverser(const geometry::LatticeData& iLatticeData) :
      mLatticeData(iLatticeData)

    {
    }

    SiteTraverser::~SiteTraverser()
    {
    }

    site_t SiteTraverser::GetXCount() const
    {
      return GetBlockSize();
    }

    site_t SiteTraverser::GetYCount() const
    {
      return GetBlockSize();
    }

    site_t SiteTraverser::GetZCount() const
    {
      return GetBlockSize();
    }

    site_t SiteTraverser::GetBlockSize() const
    {
      return mLatticeData.GetBlockSize();
    }
  }
}
