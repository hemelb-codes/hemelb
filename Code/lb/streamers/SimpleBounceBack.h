#ifndef HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H
#define HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H

#include "lb/kernels/BaseKernel.h"
#include "lb/streamers/BaseStreamer.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionType>
      class SimpleBounceBack : public BaseStreamer<SimpleBounceBack<CollisionType> >
      {
        private:
          CollisionType collider;

        public:
          SimpleBounceBack(kernels::InitParams& initParams) :
              collider(initParams)
          {

          }

          template<bool tDoRayTracing>
          inline void DoStreamAndCollide(const site_t iFirstIndex,
                                         const site_t iSiteCount,
                                         const LbmParameters* iLbmParams,
                                         geometry::LatticeData* bLatDat,
                                         hemelb::vis::Control *iControl)
          {
            for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
            {
              distribn_t *fOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);

              kernels::HydroVars<typename CollisionType::CKernel> hydroVars(fOld);

              collider.CalculatePreCollision(hydroVars, lIndex - iFirstIndex);

              collider.Collide(iLbmParams, hydroVars);

              for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
              {
                // The actual bounce-back lines, including streaming and collision. Basically swap
                // the non-equilibrium components of f in each of the opposing pairs of directions.
                site_t lStreamTo = (bLatDat->HasBoundary(lIndex, ii)) ?
                  (lIndex * D3Q15::NUMVECTORS) + D3Q15::INVERSEDIRECTIONS[ii] :
                  bLatDat->GetStreamedIndex(lIndex, ii);

                // Remember, oFNeq currently hold the equilibrium distribution. We
                // simultaneously use this and correct it, here.
                * (bLatDat->GetFNew(lStreamTo)) = hydroVars.GetFPostCollision()[ii];
              }

              //TODO: Necessary to specify sub-class?
              BaseStreamer<SimpleBounceBack>::template UpdateMinsAndMaxes<tDoRayTracing>(hydroVars.v_x,
                                                                                         hydroVars.v_y,
                                                                                         hydroVars.v_z,
                                                                                         lIndex,
                                                                                         hydroVars.GetFNeq().f,
                                                                                         hydroVars.density,
                                                                                         bLatDat,
                                                                                         iLbmParams,
                                                                                         iControl);
            }
          }

          template<bool tDoRayTracing>
          inline void DoPostStep(const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 hemelb::vis::Control *iControl)
          {

          }

          inline void DoReset(kernels::InitParams* init)
          {
            collider.Reset(init);
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_SIMPLEBOUNCEBACK_H */
