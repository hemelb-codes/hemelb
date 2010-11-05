#include "lb/collisions/ImplZeroVelocityBoundaryDensity.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      ImplZeroVelocityBoundaryDensity::ImplZeroVelocityBoundaryDensity(double* iOutletDensityArray)
      {
        mBoundaryDensityArray = iOutletDensityArray;
      }

      void ImplZeroVelocityBoundaryDensity::DoCollisions(const bool iDoRayTracing,
                                                         const double iOmega,
                                                         const int iFirstIndex,
                                                         const int iSiteCount,
                                                         MinsAndMaxes* bMinimaAndMaxima,
                                                         LocalLatticeData &bLocalLatDat,
                                                         const double iStressType,
                                                         const double iStressParam,
                                                         hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iOmega, iFirstIndex, iSiteCount,
                                      bMinimaAndMaxima, bLocalLatDat,
                                      iStressType, iStressParam, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iOmega, iFirstIndex, iSiteCount,
                                       bMinimaAndMaxima, bLocalLatDat,
                                       iStressType, iStressParam, iControl);
        }
      }

      // Collision + streaming for fluid lattice sites and adjacent to the outlet and the wall.
      template<bool tDoRayTracing>
      void ImplZeroVelocityBoundaryDensity::DoCollisionsInternal(const double iOmega,
                                                                 const int iFirstIndex,
                                                                 const int iSiteCount,
                                                                 MinsAndMaxes* bMinimaAndMaxima,
                                                                 LocalLatticeData &bLocalLatDat,
                                                                 const double iStressType,
                                                                 const double iStressParam,
                                                                 hemelb::vis::Control *iControl)
      {
        for (int iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          double *lFOld = &bLocalLatDat.FOld[iIndex * D3Q15::NUMVECTORS];
          double lFNeq[D3Q15::NUMVECTORS];
          double lDensity;

          lDensity = mBoundaryDensityArray[bLocalLatDat.GetBoundaryId(iIndex)];

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in FNeq
          D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            bLocalLatDat.FNew[bLocalLatDat.GetStreamedIndex(iIndex, ii)]
                = lFOld[ii];
            lFNeq[ii] -= lFOld[ii];
          }

          Collision::UpdateMinsAndMaxes<tDoRayTracing>(0.0, 0.0, 0.0, iIndex,
                                                       lFNeq, lDensity,
                                                       bMinimaAndMaxima,
                                                       bLocalLatDat,
                                                       iStressType,
                                                       iStressParam, iControl);
        }
      }
    }
  }
}
