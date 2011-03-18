#include "lb/collisions/ImplFInterpolation.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplFInterpolation::DoCollisions(const bool iDoRayTracing,
                                            const int iFirstIndex,
                                            const int iSiteCount,
                                            const LbmParameters &iLbmParams,
                                            MinsAndMaxes &bMinimaAndMaxima,
                                            geometry::LocalLatticeData &bLocalLatDat,
                                            hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                      bLocalLatDat, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                       bLocalLatDat, iControl);
        }
      }

      void ImplFInterpolation::PostStep(const bool iDoRayTracing,
                                        const int iFirstIndex,
                                        const int iSiteCount,
                                        const LbmParameters &iLbmParams,
                                        MinsAndMaxes &bMinimaAndMaxima,
                                        geometry::LocalLatticeData &bLocalLatDat,
                                        hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          PostStepInternal<true> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                  bLocalLatDat, iControl);
        }
        else
        {
          PostStepInternal<false> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                   bLocalLatDat, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplFInterpolation::DoCollisionsInternal(const int iFirstIndex,
                                                    const int iSiteCount,
                                                    const LbmParameters &iLbmParams,
                                                    MinsAndMaxes &bMinimaAndMaxima,
                                                    geometry::LocalLatticeData &bLocalLatDat,
                                                    hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *f = &bLocalLatDat.FOld[lIndex * D3Q15::NUMVECTORS];
          double density, v_x, v_y, v_z, f_neq[15];
          // Temporarily store f_eq in f_neq. Rectified later.
          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            bLocalLatDat.FNew[bLocalLatDat.GetStreamedIndex(lIndex, ii)] = f[ii]
                += iLbmParams.Omega * (f_neq[ii] = f[ii] - f_neq[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (v_x, v_y, v_z, lIndex, f_neq, density,
                                             bMinimaAndMaxima, bLocalLatDat, iLbmParams, iControl);
        }
      }

      //TODO: Does this change velocity / density / stress? Need to update mins and maxes if so.

      template<bool tDoRayTracing>
      void ImplFInterpolation::PostStepInternal(const int iFirstIndex,
                                                const int iSiteCount,
                                                const LbmParameters &iLbmParams,
                                                MinsAndMaxes &bMinimaAndMaxima,
                                                geometry::LocalLatticeData &bLocalLatDat,
                                                hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          // Handle odd indices, then evens - it's slightly easier to take the odd
          // and even cases separately.
          for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
          {
            if (bLocalLatDat.HasBoundary(lIndex, l))
            {
              double twoQ = 2.0 * bLocalLatDat.GetCutDistance(lIndex, l);
              if (twoQ < 1.0)
              {
                bLocalLatDat.FNew[lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]]
                    = bLocalLatDat.FNew[lIndex * D3Q15::NUMVECTORS + l] + twoQ
                        * (bLocalLatDat.FOld[lIndex * D3Q15::NUMVECTORS + l]
                            - bLocalLatDat.FNew[lIndex * D3Q15::NUMVECTORS + l]);
              }
              else
              {
                bLocalLatDat.FNew[lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]]
                    = bLocalLatDat.FOld[lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]]
                        + (1. / twoQ) * (bLocalLatDat.FOld[lIndex * D3Q15::NUMVECTORS + l]
                            - bLocalLatDat.FOld[lIndex * 15 + D3Q15::INVERSEDIRECTIONS[l]]);
              }
            }
          }
        }
      }

    }
  }
}
