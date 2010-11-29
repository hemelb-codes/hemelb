#ifndef HEMELB_LB_LOCALLATTICEDATA_H
#define HEMELB_LB_LOCALLATTICEDATA_H

#include "constants.h"
#include "D3Q15.h"
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

    class LocalLatticeData
    {
      public:
        LocalLatticeData(int iLocalFluidSites, int iSharedFluidSites)
        {
          LocalFluidSites = iLocalFluidSites;

          // Allocate f_old and f_new according to the number of sites on the process.  The extra site
          // is there for when we would stream into a solid site during the simulation, which avoids
          // an if condition at every timestep at every boundary site.  We also allocate space for the
          // shared distribution functions.  We need twice as much space when we check the convergence
          // and the extra distribution functions are
          FOld = new double[LocalFluidSites * D3Q15::NUMVECTORS + 1
              + iSharedFluidSites];
          FNew = new double[LocalFluidSites * D3Q15::NUMVECTORS + 1
              + iSharedFluidSites];

          // f_id is allocated so we know which sites to get information from.
          mFNeighbours = new int[LocalFluidSites * D3Q15::NUMVECTORS];

          mSiteData = new unsigned int[LocalFluidSites];
          mWallNormalAtSite = new double[LocalFluidSites * 3];
          mDistanceToWall = new double[LocalFluidSites
              * (D3Q15::NUMVECTORS - 1)];
        }

        ~LocalLatticeData()
        {
          delete[] FOld;
          delete[] FNew;
          delete[] mSiteData;
          delete[] mFNeighbours;
          delete[] mDistanceToWall;
          delete[] mWallNormalAtSite;
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

        const double *GetNormalToWall(int iSiteIndex) const
        {
          return &mWallNormalAtSite[iSiteIndex * 3];
        }

        SiteType GetSiteType(int iSiteIndex) const
        {
          return (SiteType) (mSiteData[iSiteIndex] & SITE_TYPE_MASK);
        }

        int GetLocalFluidSiteCount() const
        {
          return LocalFluidSites;
        }

        void SetNeighbourLocation(int iSiteIndex, int iDirection, int iValue)
        {
          mFNeighbours[iSiteIndex * D3Q15::NUMVECTORS + iDirection] = iValue;
        }

        void SetWallNormal(int iSiteIndex, const double iNormal[3])
        {
          for (int ii = 0; ii < 3; ii++)
            mWallNormalAtSite[iSiteIndex * 3 + ii] = iNormal[ii];
        }

        void SetDistanceToWall(int iSiteIndex,
                               const double iCutDistance[D3Q15::NUMVECTORS - 1])
        {
          for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
            mDistanceToWall[iSiteIndex * (D3Q15::NUMVECTORS - 1) + l]
                = iCutDistance[l];
        }

      public:
        double *FOld;
        double *FNew;

        // TODO sadly this has to be public, due to some budgetry in the way we determine site type.
        // SiteType || FluidSite and SiteType && FluidSite have different significances...
        unsigned int *mSiteData;

      private:
        int LocalFluidSites;
        int *mFNeighbours;
        double *mDistanceToWall;
        double *mWallNormalAtSite;
    };
  }
}

#endif /* HEMELB_LB_LOCALLATTICEDATA_H */
