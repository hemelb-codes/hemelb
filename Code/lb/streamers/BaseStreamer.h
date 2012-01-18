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
          inline void StreamAndCollide(const site_t iFirstIndex,
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
          inline void PostStep(const site_t iFirstIndex,
                               const site_t iSiteCount,
                               const LbmParameters* iLbmParams,
                               geometry::LatticeData* bLatDat,
                               hemelb::vis::Control *iControl)
          {
            // The template parameter is required because we're using the CRTP to call a
            // metaprogrammed method of the implementation class.
            static_cast<StreamerImpl*>(this)->template DoPostStep<tDoRayTracing>(iFirstIndex,
                                                                                 iSiteCount,
                                                                                 iLbmParams,
                                                                                 bLatDat,
                                                                                 iControl);
          }

          inline void Reset(kernels::InitParams* init)
          {
            static_cast<StreamerImpl*>(this)->DoReset(init);
          }

        protected:
          template<bool tDoRayTracing>
          inline static void UpdateMinsAndMaxes(distribn_t iVx,
                                                distribn_t iVy,
                                                distribn_t iVz,
                                                const geometry::Site& site,
                                                const distribn_t* f_neq,
                                                const distribn_t iDensity,
                                                const LbmParameters* iLbmParams,
                                                hemelb::vis::Control *iControl)
          {
            if (tDoRayTracing)
            {
              distribn_t rtStress;

              if (iLbmParams->StressType == ShearStress)
              {
                if (!site.IsEdge())
                {
                  rtStress = NO_VALUE;
                }
                else
                {
                  D3Q15::CalculateShearStress(iDensity,
                                              f_neq,
                                              site.GetWallNormal(),
                                              rtStress,
                                              iLbmParams->GetStressParameter());
                }
              }
              else
              {
                D3Q15::CalculateVonMisesStress(f_neq, rtStress, iLbmParams->GetStressParameter());
              }

              // TODO: It'd be nice if the /iDensity were unnecessary.
              distribn_t lVelocity = sqrt(iVx * iVx + iVy * iVy + iVz * iVz) / iDensity;
              iControl->RegisterSite(site.GetIndex(), iDensity, lVelocity, rtStress);
            }
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BASESTREAMER_H */
