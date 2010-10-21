#include "lb/collisions/ImplFInterpolation.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplFInterpolation::DoCollisions(const bool iDoRayTracing,
                                            const double iOmega,
                                            double iFOldAll[],
                                            double iFNewAll[],
                                            const int iFIdAll[],
                                            const int iFirstIndex,
                                            const int iSiteCount,
                                            MinsAndMaxes* bMinimaAndMaxima,
                                            const Net* net,
                                            const double iStressType,
                                            const double iStressParam,
                                            hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                      iFirstIndex, iSiteCount,
                                      bMinimaAndMaxima, net, iStressType,
                                      iStressParam, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                       iFirstIndex, iSiteCount,
                                       bMinimaAndMaxima, net, iStressType,
                                       iStressParam, iControl);
        }
      }

      void ImplFInterpolation::PostStep(const bool iDoRayTracing,
                                        const double iOmega,
                                        double iFOldAll[],
                                        double iFNewAll[],
                                        const int iFIdAll[],
                                        const int iFirstIndex,
                                        const int iSiteCount,
                                        MinsAndMaxes* bMinimaAndMaxima,
                                        const Net* net,
                                        const double iStressType,
                                        const double iStressParam,
                                        hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          PostStepInternal<true> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                  iFirstIndex, iSiteCount, bMinimaAndMaxima,
                                  net, iStressType, iStressParam, iControl);
        }
        else
        {
          PostStepInternal<false> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                   iFirstIndex, iSiteCount, bMinimaAndMaxima,
                                   net, iStressType, iStressParam, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplFInterpolation::DoCollisionsInternal(const double iOmega,
                                                    double iFOldAll[],
                                                    double iFNewAll[],
                                                    const int iFIdAll[],
                                                    const int iFirstIndex,
                                                    const int iSiteCount,
                                                    MinsAndMaxes* bMinimaAndMaxima,
                                                    const Net* net,
                                                    const double iStressType,
                                                    const double iStressParam,
                                                    hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *f = &iFOldAll[lIndex * D3Q15::NUMVECTORS];
          double density, v_x, v_y, v_z, f_neq[15];
          // Temporarily store f_eq in f_neq. Rectified later.
          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            iFNewAll[iFIdAll[lIndex * D3Q15::NUMVECTORS + ii]] = f[ii]
                += iOmega * (f_neq[ii] = f[ii] - f_neq[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (v_x, v_y, v_z, lIndex, f_neq,
                                             density, bMinimaAndMaxima, net,
                                             iStressType, iStressParam, iControl);
        }
      }

      //TODO: Does this change velocity / density / stress? Need to update mins and maxes if so.

      template<bool tDoRayTracing>
      void ImplFInterpolation::PostStepInternal(const double iOmega,
                                                double iFOldAll[],
                                                double iFNewAll[],
                                                const int iFIdAll[],
                                                const int iFirstIndex,
                                                const int iSiteCount,
                                                MinsAndMaxes* bMinimaAndMaxima,
                                                const Net* net,
                                                const double iStressType,
                                                const double iStressParam,
                                                hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          // Handle odd indices, then evens - it's slightly easier to take the odd
          // and even cases separately.
          for (int l = 1; l < D3Q15::NUMVECTORS; l++)
          {
            if (net->HasBoundary(lIndex, l))
            {
              double twoQ = 2.0 * net->GetCutDistance(lIndex, l);
              if (twoQ < 1.0)
              {
                f_new[lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]]
                    = f_new[lIndex * D3Q15::NUMVECTORS + l] + twoQ
                        * (f_old[lIndex * D3Q15::NUMVECTORS + l] - f_new[lIndex
                            * D3Q15::NUMVECTORS + l]);
              }
              else
              {
                f_new[lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]]
                    = f_old[lIndex * D3Q15::NUMVECTORS
                        + D3Q15::INVERSEDIRECTIONS[l]] + (1. / twoQ)
                        * (f_old[lIndex * D3Q15::NUMVECTORS + l] - f_old[lIndex
                            * 15 + D3Q15::INVERSEDIRECTIONS[l]]);
              }
            }
          }
        }
      }

    }
  }
}
