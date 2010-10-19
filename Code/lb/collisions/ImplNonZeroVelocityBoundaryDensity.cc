#include "lb/collisions/ImplNonZeroVelocityBoundaryDensity.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      ImplNonZeroVelocityBoundaryDensity::ImplNonZeroVelocityBoundaryDensity(double* iOutletDensityArray)
      {
        mBoundaryDensityArray = iOutletDensityArray;
      }

      void ImplNonZeroVelocityBoundaryDensity::DoCollisions(const bool iDoRayTracing,
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

      template<bool tDoRayTracing>
      void ImplNonZeroVelocityBoundaryDensity::DoCollisionsInternal(const double iOmega,
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
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *lFOld = &iFOldAll[lIndex * D3Q15::NUMVECTORS];
          double lFNeq[15];
          double lVx, lVy, lVz, lDummyDensity, lDensity;

          lDensity = mBoundaryDensityArray[net->GetBoundaryId(lIndex)];

          D3Q15::CalculateDensityAndVelocity(lFOld, lDummyDensity, lVx, lVy,
                                             lVz);

          for(unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq (rectified later).
          D3Q15::CalculateFeq(lDensity, lVx, lVy, lVz, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            iFNewAll[iFIdAll[lIndex * D3Q15::NUMVECTORS + ii]] = lFOld[ii];
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] -= lFOld[ii];
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx, lVy, lVz, lIndex, lFNeq,
                                             lDensity, bMinimaAndMaxima, net,
                                             iStressType, iStressParam);
        }
      }

    }
  }
}
