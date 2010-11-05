#ifndef HEMELB_LB_LOCALLATTICEDATA_H
#define HEMELB_LB_LOCALLATTICEDATA_H

#include "debug/Debugger.h"

namespace hemelb
{
  namespace lb
  {
    enum SiteType
    {
      SOLID_TYPE = 0U,
      FLUID_TYPE = 1U,
      INLET_TYPE = 2U,
      OUTLET_TYPE = 3U
    };

    struct LocalLatticeData
    {
      public:
        LocalLatticeData()
        {
          FOld = NULL;
          FNew = NULL;
          mSiteData = NULL;
          mFNeighbours = NULL;
          mDistanceToWall = NULL;
          mWallNormalAtSite = NULL;
        }

        ~LocalLatticeData()
        {
          hemelb::debug::Debugger::Get()->BreakHere();

          if (FOld != NULL)
          {
            delete[] FOld;
          }
          if (FNew != NULL)
          {
            delete[] FNew;
          }
          if (mSiteData != NULL)
          {
            delete[] mSiteData;
          }
          if (mFNeighbours != NULL)
          {
            delete[] mFNeighbours;
          }
          if (mDistanceToWall != NULL)
          {
            delete[] mDistanceToWall;
          }
          if (mWallNormalAtSite != NULL)
          {
            delete[] mWallNormalAtSite;
          }
        }

        int GetStreamedIndex(int iSiteIndex, int iDirectionIndex) const
        {
          return mFNeighbours[iSiteIndex * D3Q15::NUMVECTORS + iDirectionIndex];
        }

        double GetCutDistance(int iSiteIndex, int iDirection) const
        {
          return mDistanceToWall[iSiteIndex * (D3Q15::NUMVECTORS - 1)
              + iDirection - 1];
        }

        bool HasBoundary(int iSiteIndex, int iDirection) const
        {
          unsigned int lBoundaryConfig = (mSiteData[iSiteIndex]
              & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT;
          return (lBoundaryConfig & (1U << (iDirection - 1))) != 0;
        }

        int GetBoundaryId(int iSiteIndex) const
        {
          return (mSiteData[iSiteIndex] & BOUNDARY_ID_MASK)
              >> BOUNDARY_ID_SHIFT;
        }

        double *GetNormalToWall(int iSiteIndex) const
        {
          return &mWallNormalAtSite[iSiteIndex * 3];
        }

        SiteType GetSiteType(int iSiteIndex)
        {
          return (SiteType) (mSiteData[iSiteIndex] & SITE_TYPE_MASK);
        }

        double *FOld;
        double *FNew;

        friend class Net;

        //private:
        // TODO In time, these need to be made private. At the moment, that's impractical because they're initialised
        // by the Net class.
        unsigned int *mSiteData;
        int *mFNeighbours;
        double *mDistanceToWall;
        double *mWallNormalAtSite;
    };
  }
}

#endif /* HEMELB_LB_LOCALLATTICEDATA_H */
