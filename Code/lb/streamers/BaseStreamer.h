#ifndef HEMELB_LB_STREAMERS_BASESTREAMER_H
#define HEMELB_LB_STREAMERS_BASESTREAMER_H

#include <math.h>

#include "geometry/LatticeData.h"
#include "vis/Control.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      /**
       * BaseStreamer: inheritable base class for the streaming operator. The public interface
       * here defines the complete interface usable by external code.
       *  - Constructor(InitParams&)
       *  - <bool tDoRayTracing> StreamAndCollide(const site_t, const site_t, const LbmParameters*,
       *      geometry::LatticeData*, hemelb::vis::Control*)
       *  - <bool tDoRayTracing> PostStep(const site_t, const site_t, const LbmParameters*,
       *      geometry::LatticeData*, hemelb::vis::Control*)
       *  - Reset(kernels::InitParams* init)
       *
       * The following must be implemented by concrete streamers (which derive from this class
       * using the CRTP).
       *  - Constructor(InitParams&)
       *  - <bool tDoRayTracing> DoStreamAndCollide(const site_t, const site_t, const LbmParameters*,
       *      geometry::LatticeData*, hemelb::vis::Control*)
       *  - <bool tDoRayTracing> DoPostStep(const site_t, const site_t, const LbmParameters*,
       *      geometry::LatticeData*, hemelb::vis::Control*)
       *  - DoReset(kernels::InitParams* init)
       */
      template<typename StreamerImpl>
      class BaseStreamer
      {
        public:
          template<bool tDoRayTracing>
          void StreamAndCollide(const site_t iFirstIndex,
                                const site_t iSiteCount,
                                const LbmParameters* iLbmParams,
                                geometry::LatticeData* bLatDat,
                                hemelb::vis::Control *iControl)
          {
            static_cast<StreamerImpl*>(this)->template DoStreamAndCollide<tDoRayTracing>(iFirstIndex,
                                                                                         iSiteCount,
                                                                                         iLbmParams,
                                                                                         bLatDat,
                                                                                         iControl);
          }

          template<bool tDoRayTracing>
          void PostStep(const site_t iFirstIndex,
                        const site_t iSiteCount,
                        const LbmParameters* iLbmParams,
                        geometry::LatticeData* bLatDat,
                        hemelb::vis::Control *iControl)
          {
            static_cast<StreamerImpl*>(this)->DoPostStep<tDoRayTracing>(iFirstIndex,
                                                                        iSiteCount,
                                                                        iLbmParams,
                                                                        bLatDat,
                                                                        iControl);
          }

          void Reset(kernels::InitParams* init)
          {
            static_cast<StreamerImpl*>(this)->DoReset(init);
          }

        protected:
          template<bool tDoRayTracing>
          static void UpdateMinsAndMaxes(distribn_t iVx,
                                         distribn_t iVy,
                                         distribn_t iVz,
                                         const site_t iSiteIndex,
                                         const distribn_t* f_neq,
                                         const distribn_t iDensity,
                                         const geometry::LatticeData* iLatDat,
                                         const LbmParameters* iLbmParams,
                                         hemelb::vis::Control *iControl)
          {
            if (tDoRayTracing)
            {
              distribn_t rtStress;

              if (iLbmParams->StressType == ShearStress)
              {
                if (iLatDat->GetNormalToWall(iSiteIndex)[0] > NO_VALUE)
                {
                  rtStress = NO_VALUE;
                }
                else
                {
                  D3Q15::CalculateShearStress(iDensity,
                                              f_neq,
                                              iLatDat->GetNormalToWall(iSiteIndex),
                                              rtStress,
                                              iLbmParams->StressParameter());
                }
              }
              else
              {
                D3Q15::CalculateVonMisesStress(f_neq, rtStress, iLbmParams->StressParameter());
              }

              // TODO: It'd be nice if the /iDensity were unnecessary.
              distribn_t lVelocity = sqrt(iVx * iVx + iVy * iVy + iVz * iVz) / iDensity;
              iControl->RegisterSite(iSiteIndex, iDensity, lVelocity, rtStress);
            }
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BASESTREAMER_H */
