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
                                                         double iFOldAll[],
                                                         double iFNewAll[],
                                                         const int iFIdAll[],
                                                         const int iFirstIndex,
                                                         const int iSiteCount,
                                                         MinsAndMaxes* bMinimaAndMaxima,
                                                         const Net* net,
                                                         const double iStressType,
                                                         const double iStressParam)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                      iFirstIndex, iSiteCount,
                                      bMinimaAndMaxima, net, iStressType,
                                      iStressParam);
        }
        else
        {
          DoCollisionsInternal<false> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                       iFirstIndex, iSiteCount,
                                       bMinimaAndMaxima, net, iStressType,
                                       iStressParam);
        }
      }

      // Collision + streaming for fluid lattice sites and adjacent to the outlet and the wall.
      template<bool tDoRayTracing>
      void ImplZeroVelocityBoundaryDensity::DoCollisionsInternal(const double iOmega,
                                                                 double iFOldAll[],
                                                                 double iFNewAll[],
                                                                 const int iFIdAll[],
                                                                 const int iFirstIndex,
                                                                 const int iSiteCount,
                                                                 MinsAndMaxes* bMinimaAndMaxima,
                                                                 const Net* net,
                                                                 const double iStressType,
                                                                 const double iStressParam)
      {
        for (int iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          double *lFOld = &iFOldAll[iIndex * D3Q15::NUMVECTORS];
          double lFNeq[15];
          double lDensity;

          lDensity = mBoundaryDensityArray[net->GetBoundaryId(iIndex)];

          for(int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in FNeq
          D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            iFNewAll[iFIdAll[iIndex * D3Q15::NUMVECTORS + ii]] = lFOld[ii];
            lFNeq[ii]-= lFOld[ii];
          }

          Collision::UpdateMinsAndMaxes<tDoRayTracing> (0.0, 0.0, 0.0, iIndex, lFNeq,
                                             lDensity, bMinimaAndMaxima, net,
                                             iStressType, iStressParam);
        }
      }
    }
  }
}
