#include "lb/collisions/ImplSimpleBounceBack.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplSimpleBounceBack::DoCollisions(const bool iDoRayTracing,
                                              const int iFirstIndex,
                                              const int iSiteCount,
                                              const LbmParameters &iLbmParams,
                                              MinsAndMaxes* bMinimaAndMaxima,
                                              LocalLatticeData &bLocalLatDat,
                                              hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iFirstIndex, iSiteCount, iLbmParams,
                                      bMinimaAndMaxima, bLocalLatDat, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iFirstIndex, iSiteCount, iLbmParams,
                                       bMinimaAndMaxima, bLocalLatDat, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplSimpleBounceBack::DoCollisionsInternal(const int iFirstIndex,
                                                      const int iSiteCount,
                                                      const LbmParameters &iLbmParams,
                                                      MinsAndMaxes* bMinimaAndMaxima,
                                                      LocalLatticeData &bLocalLatDat,
                                                      hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *lFOld = &bLocalLatDat.FOld[lIndex * D3Q15::NUMVECTORS];
          double lFNeq[D3Q15::NUMVECTORS];
          double lVx, lVy, lVz, lDensity;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          D3Q15::CalculateDensityVelocityFEq(lFNeq, lDensity, lVx, lVy, lVz,
                                             lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            // The actual bounce-back lines, including streaming and collision. Basically swap the non-equilibrium components of f in each of the opposing pairs of directions.
            int lStreamTo = (bLocalLatDat.HasBoundary(lIndex, ii))
              ? ( (lIndex * D3Q15::NUMVECTORS) + D3Q15::INVERSEDIRECTIONS[ii])
              : bLocalLatDat.GetStreamedIndex(lIndex, ii);

            // Remember, oFNeq currently hold the equilibrium distribution. We
            // simultaneously use this and correct it, here.
            bLocalLatDat.FNew[lStreamTo] = lFOld[ii] += iLbmParams.Omega
                * (lFNeq[ii] -= lFOld[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx, lVy, lVz, lIndex, lFNeq,
                                             lDensity, bMinimaAndMaxima,
                                             bLocalLatDat, iLbmParams, iControl);
        }
      }

    }
  }
}
