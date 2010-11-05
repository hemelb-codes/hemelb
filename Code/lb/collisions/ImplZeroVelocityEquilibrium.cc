#include "lb/collisions/ImplZeroVelocityEquilibrium.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      void ImplZeroVelocityEquilibrium::DoCollisions(const bool iDoRayTracing,
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

      template<bool tDoRayTracing>
      void ImplZeroVelocityEquilibrium::DoCollisionsInternal(const double iOmega,
                                                             const int iFirstIndex,
                                                             const int iSiteCount,
                                                             MinsAndMaxes* bMinimaAndMaxima,
                                                             LocalLatticeData &bLocalLatDat,
                                                             const double iStressType,
                                                             const double iStressParam,
                                                             hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *lFOld = &bLocalLatDat.FOld[lIndex * D3Q15::NUMVECTORS];
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
            bLocalLatDat.FNew[bLocalLatDat.GetStreamedIndex(lIndex, ii)]
                = lFOld[ii];
            lFNeq[ii] -= lFOld[ii];
          }

          Collision::UpdateMinsAndMaxes<tDoRayTracing>(0.0, 0.0, 0.0, lIndex,
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
