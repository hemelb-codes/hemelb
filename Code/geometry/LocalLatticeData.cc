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

    void LatticeData::LocalLatticeData::Initialise(int iLocalFluidSites)
    {
      LocalFluidSites = iLocalFluidSites;

      // f_id is allocated so we know which sites to get information from.
      mFNeighbours = new unsigned int[LocalFluidSites * D3Q15::NUMVECTORS];

      mSiteData = new unsigned int[LocalFluidSites];
      mWallNormalAtSite = new double[LocalFluidSites * 3];
      mDistanceToWall = new double[LocalFluidSites * (D3Q15::NUMVECTORS - 1)];
    }

    void LatticeData::LocalLatticeData::SetSharedSiteCount(int iSharedCount)
    {
      // Allocate f_old and f_new according to the number of sites on the process.  The extra site
      // is there for when we would stream into a solid site during the simulation, which avoids
      // an if condition at every timestep at every boundary site.  We also allocate space for the
      // shared distribution functions.  We need twice as much space when we check the convergence
      // and the extra distribution functions are
      FOld = new double[LocalFluidSites * D3Q15::NUMVECTORS + 1 + iSharedCount];
      FNew = new double[LocalFluidSites * D3Q15::NUMVECTORS + 1 + iSharedCount];
    }

    LatticeData::LocalLatticeData::~LocalLatticeData()
    {
      if (FOld != NULL)
        delete[] FOld;
      if (FNew != NULL)
        delete[] FNew;
      if (mSiteData != NULL)
        delete[] mSiteData;
      if (mFNeighbours != NULL)
        delete[] mFNeighbours;
      if (mDistanceToWall != NULL)
        delete[] mDistanceToWall;
      if (mWallNormalAtSite != NULL)
        delete[] mWallNormalAtSite;
    }

    int LatticeData::LocalLatticeData::GetStreamedIndex(int iSiteIndex, int iDirectionIndex) const
    {
      return mFNeighbours[iSiteIndex * D3Q15::NUMVECTORS + iDirectionIndex];
    }

    double LatticeData::LocalLatticeData::GetCutDistance(int iSiteIndex, int iDirection) const
    {
      return mDistanceToWall[iSiteIndex * (D3Q15::NUMVECTORS - 1) + iDirection - 1];
    }

    bool LatticeData::LocalLatticeData::HasBoundary(int iSiteIndex, int iDirection) const
    {
      unsigned int lBoundaryConfig = (mSiteData[iSiteIndex] & BOUNDARY_CONFIG_MASK)
          >> BOUNDARY_CONFIG_SHIFT;
      return (lBoundaryConfig & (1U << (iDirection - 1))) != 0;
    }

    int LatticeData::LocalLatticeData::GetBoundaryId(int iSiteIndex) const
    {
      return (mSiteData[iSiteIndex] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
    }

    const double* LatticeData::LocalLatticeData::GetNormalToWall(int iSiteIndex) const
    {
      return &mWallNormalAtSite[iSiteIndex * 3];
    }

    LatticeData::SiteType LatticeData::LocalLatticeData::GetSiteType(int iSiteIndex) const
    {
      return (SiteType) (mSiteData[iSiteIndex] & SITE_TYPE_MASK);
    }

    int LatticeData::LocalLatticeData::GetLocalFluidSiteCount() const
    {
      return LocalFluidSites;
    }

    void LatticeData::LocalLatticeData::SetNeighbourLocation(unsigned int iSiteIndex,
                                                             unsigned int iDirection,
                                                             unsigned int iValue)
    {
      mFNeighbours[iSiteIndex * D3Q15::NUMVECTORS + iDirection] = iValue;
    }

    void LatticeData::LocalLatticeData::SetWallNormal(int iSiteIndex, const double iNormal[3])
    {
      for (int ii = 0; ii < 3; ii++)
        mWallNormalAtSite[iSiteIndex * 3 + ii] = iNormal[ii];
    }

    void LatticeData::LocalLatticeData::SetDistanceToWall(int iSiteIndex,
                                                          const double iCutDistance[D3Q15::NUMVECTORS
                                                              - 1])
    {
      for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
        mDistanceToWall[iSiteIndex * (D3Q15::NUMVECTORS - 1) + l] = iCutDistance[l];
    }

  }
}
