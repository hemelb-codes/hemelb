#ifndef HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H
#define HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H

#include "lb/streamers/BaseStreamer.h"

/**
 * TODO This class hasn't been modified to fit the new kernel / collision / streamer
 * hierarchy. This needs to happen before it can be used in the LBM.
 */
namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename tCollisionOperator>
      class SimpleBounceBack : public Implementation
      {
        public:
          template<bool tDoRayTracing>
          void DoStreamAndCollide(WallCollision* mWallCollision,
                                  const site_t iFirstIndex,
                                  const site_t iSiteCount,
                                  const LbmParameters* iLbmParams,
                                  geometry::LatticeData* bLatDat,
                                  hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void DoPostStep(WallCollision* mWallCollision,
                          const site_t iFirstIndex,
                          const site_t iSiteCount,
                          const LbmParameters* iLbmParams,
                          geometry::LatticeData* bLatDat,
                          hemelb::vis::Control *iControl);

      };

      template<typename tCollisionOperator>
      template<bool tDoRayTracing>
      void SimpleBounceBack<tCollisionOperator>::DoStreamAndCollide(WallCollision* mWallCollision,
                                                                    const site_t iFirstIndex,
                                                                    const site_t iSiteCount,
                                                                    const LbmParameters* iLbmParams,
                                                                    geometry::LatticeData* bLatDat,
                                                                    hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t *lFOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t lFNeq[D3Q15::NUMVECTORS];
          distribn_t lVx, lVy, lVz, lDensity;

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          D3Q15::CalculateDensityVelocityFEq(lFNeq, lDensity, lVx, lVy, lVz, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            // The actual bounce-back lines, including streaming and collision. Basically swap the non-equilibrium components of f in each of the opposing pairs of directions.
            site_t lStreamTo = (bLatDat->HasBoundary(lIndex, ii))
            ?
              ( (lIndex * D3Q15::NUMVECTORS) + D3Q15::INVERSEDIRECTIONS[ii])
              :
              bLatDat->GetStreamedIndex(lIndex, ii);

            // Remember, oFNeq currently hold the equilibrium distribution. We
            // simultaneously use this and correct it, here.
            * (bLatDat->GetFNew(lStreamTo)) = lFOld[ii] += iLbmParams->Omega * (lFNeq[ii] -=
                lFOld[ii]);
          }

          UpdateMinsAndMaxes < tDoRayTracing
              > (lVx, lVy, lVz, lIndex, lFNeq, lDensity, bLatDat, iLbmParams, iControl);
        }
      }

      template<typename tCollisionOperator>
      template<bool tDoRayTracing>
      void SimpleBounceBack<tCollisionOperator>::DoPostStep(WallCollision* mWallCollision,
                                                            const site_t iFirstIndex,
                                                            const site_t iSiteCount,
                                                            const LbmParameters* iLbmParams,
                                                            geometry::LatticeData* bLatDat,
                                                            hemelb::vis::Control *iControl)
      {

      }

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H */
