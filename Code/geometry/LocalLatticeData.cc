#include "LatticeData.h"

namespace hemelb
{
  namespace geometry
  {

    LatticeData::LocalLatticeData::LocalLatticeData()
    {
      FOld = NULL;
      FNew = NULL;
      mSiteData = NULL;
      mFNeighbours = NULL;
      mDistanceToWall = NULL;
      mWallNormalAtSite = NULL;
    }

    LatticeData::LocalLatticeData::LocalLatticeData(site_t iLocalFluidSites)
    {
      LocalFluidSites = iLocalFluidSites;

      // f_id is allocated so we know which sites to get information from.
      mFNeighbours = new site_t[LocalFluidSites * D3Q15::NUMVECTORS];

      mSiteData = new unsigned int[LocalFluidSites];
      mWallNormalAtSite = new double[LocalFluidSites * 3];
      mDistanceToWall = new double[LocalFluidSites * (D3Q15::NUMVECTORS - 1)];

      FOld = NULL;
      FNew = NULL;
    }

    void LatticeData::LocalLatticeData::SetSharedSiteCount(site_t iSharedCount)
    {
      // Allocate f_old and f_new according to the number of sites on the process.  The extra site
      // is there for when we would stream into a solid site during the simulation, which avoids
      // an if condition at every timestep at every boundary site.  We also allocate space for the
      // shared distribution functions.  We need twice as much space when we check the convergence
      // and the extra distribution functions are
      FOld = new distribn_t[LocalFluidSites * D3Q15::NUMVECTORS + 1 + iSharedCount];
      FNew = new distribn_t[LocalFluidSites * D3Q15::NUMVECTORS + 1 + iSharedCount];
    }

    LatticeData::LocalLatticeData::~LocalLatticeData()
    {
      delete[] FOld;
      delete[] FNew;
      delete[] mSiteData;
      delete[] mFNeighbours;
      delete[] mDistanceToWall;
      delete[] mWallNormalAtSite;
    }

    site_t LatticeData::LocalLatticeData::GetStreamedIndex(site_t iSiteIndex,
                                                           unsigned int iDirectionIndex) const
    {
      return mFNeighbours[iSiteIndex * D3Q15::NUMVECTORS + iDirectionIndex];
    }

    double LatticeData::LocalLatticeData::GetCutDistance(site_t iSiteIndex, int iDirection) const
    {
      return mDistanceToWall[iSiteIndex * (D3Q15::NUMVECTORS - 1) + iDirection - 1];
    }

    bool LatticeData::LocalLatticeData::HasBoundary(const site_t siteIndex, const int direction) const
    {
      /*
       * Have a boundary in the specified direction if the corresponding bit of the boundary config is set.
       * Direction 0 points to the current site, hence there can be no boundary.
       */
      if (direction > 0)
      {
        const unsigned int directionMask = 1U << (direction - 1);
        const unsigned int boundaryBits = (mSiteData[siteIndex] & BOUNDARY_CONFIG_MASK);
        const unsigned int shiftedBoundaryBits = boundaryBits >> BOUNDARY_CONFIG_SHIFT;
        return (shiftedBoundaryBits & directionMask) != 0;
      }
      else
      {
        return false;
      }
    }

    int LatticeData::LocalLatticeData::GetBoundaryId(site_t iSiteIndex) const
    {
      return (mSiteData[iSiteIndex] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
    }

    const double* LatticeData::LocalLatticeData::GetNormalToWall(site_t iSiteIndex) const
    {
      return &mWallNormalAtSite[iSiteIndex * 3];
    }

    LatticeData::SiteType LatticeData::LocalLatticeData::GetSiteType(site_t iSiteIndex) const
    {
      return (SiteType) (mSiteData[iSiteIndex] & SITE_TYPE_MASK);
    }

    site_t LatticeData::LocalLatticeData::GetLocalFluidSiteCount() const
    {
      return LocalFluidSites;
    }

    void LatticeData::LocalLatticeData::SetNeighbourLocation(site_t iSiteIndex,
                                                             unsigned int iDirection,
                                                             site_t iValue)
    {
      mFNeighbours[iSiteIndex * D3Q15::NUMVECTORS + iDirection] = iValue;
    }

    void LatticeData::LocalLatticeData::SetWallNormal(site_t iSiteIndex, const double iNormal[3])
    {
      for (int ii = 0; ii < 3; ii++)
        mWallNormalAtSite[iSiteIndex * 3 + ii] = iNormal[ii];
    }

    void LatticeData::LocalLatticeData::SetDistanceToWall(site_t iSiteIndex,
                                                          const double iCutDistance[D3Q15::NUMVECTORS
                                                              - 1])
    {
      for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
        mDistanceToWall[iSiteIndex * (D3Q15::NUMVECTORS - 1) + l] = iCutDistance[l];
    }

  }
}
