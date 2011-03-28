#include "lb/collisions/ImplZeroVelocityEquilibrium.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      void ImplZeroVelocityEquilibrium::DoCollisions(const bool iDoRayTracing,
                                                     const int iFirstIndex,
                                                     const int iSiteCount,
                                                     const LbmParameters &iLbmParams,
                                                     MinsAndMaxes &bMinimaAndMaxima,
                                                     geometry::LatticeData &bLatDat,
                                                     hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                      bLatDat, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                       bLatDat, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplZeroVelocityEquilibrium::DoCollisionsInternal(const int iFirstIndex,
                                                             const int iSiteCount,
                                                             const LbmParameters &iLbmParams,
                                                             MinsAndMaxes &bMinimaAndMaxima,
                                                             geometry::LatticeData &bLatDat,
                                                             hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *lFOld = bLatDat.GetFOld(lIndex * D3Q15::NUMVECTORS);
          double lFNeq[D3Q15::NUMVECTORS];
          double lDensity;

          lDensity = 0.0;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lDensity += lFOld[ii];
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq
          D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat.GetFNew(bLatDat.GetStreamedIndex(lIndex, ii))) = lFOld[ii];
            lFNeq[ii] -= lFOld[ii];
          }

          Collision::UpdateMinsAndMaxes<tDoRayTracing>(0.0, 0.0, 0.0, lIndex, lFNeq, lDensity,
                                                       bMinimaAndMaxima, bLatDat, iLbmParams,
                                                       iControl);
        }
      }
    }
  }
}
