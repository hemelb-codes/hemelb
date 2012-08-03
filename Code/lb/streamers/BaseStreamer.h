// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_STREAMERS_BASESTREAMER_H
#define HEMELB_LB_STREAMERS_BASESTREAMER_H

#include <cmath>

#include "geometry/LatticeData.h"
#include "vis/Control.h"
#include "lb/LbmParameters.h"
#include "lb/MacroscopicPropertyCache.h"

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
       *  - typedef for CollisionType, the type of the collider operation.
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
          inline void StreamAndCollide(const site_t firstIndex,
                                       const site_t siteCount,
                                       const LbmParameters* lbmParams,
                                       geometry::LatticeData* latDat,
                                       lb::MacroscopicPropertyCache& propertyCache)
          {
            static_cast<StreamerImpl*>(this)->template DoStreamAndCollide<tDoRayTracing>(firstIndex,
                                                                                         siteCount,
                                                                                         lbmParams,
                                                                                         latDat,
                                                                                         propertyCache);
          }

          template<bool tDoRayTracing>
          inline void PostStep(const site_t firstIndex,
                               const site_t siteCount,
                               const LbmParameters* lbmParams,
                               geometry::LatticeData* latDat,
                               lb::MacroscopicPropertyCache& propertyCache)
          {
            // The template parameter is required because we're using the CRTP to call a
            // metaprogrammed method of the implementation class.
            static_cast<StreamerImpl*>(this)->template DoPostStep<tDoRayTracing>(firstIndex,
                                                                                 siteCount,
                                                                                 lbmParams,
                                                                                 latDat,
                                                                                 propertyCache);
          }

          inline void Reset(kernels::InitParams* init)
          {
            static_cast<StreamerImpl*>(this)->DoReset(init);
          }

        protected:
          template<bool tDoRayTracing>
          inline static void UpdateMinsAndMaxes(distribn_t velocity_x,
                                                distribn_t velocity_y,
                                                distribn_t velocity_z,
                                                const geometry::Site& site,
                                                const distribn_t* f_neq,
                                                const distribn_t density,
                                                const distribn_t tau,
                                                const LbmParameters* lbmParams,
                                                lb::MacroscopicPropertyCache& propertyCache)
          {
            if (propertyCache.densityCache.RequiresRefresh())
            {
              propertyCache.densityCache.Put(site.GetIndex(), density);
            }

            if (propertyCache.velocityCache.RequiresRefresh())
            {
              propertyCache.velocityCache.Put(site.GetIndex(),
                                              util::Vector3D<distribn_t>(velocity_x, velocity_y, velocity_z));
            }

            if (propertyCache.shearStressCache.RequiresRefresh())
            {
              distribn_t stress;

              if (!site.IsEdge())
              {
                stress = NO_VALUE;
              }
              else
              {
                StreamerImpl::CollisionType::CKernel::LatticeType::CalculateShearStress(density,
                                                                                        f_neq,
                                                                                        site.GetWallNormal(),
                                                                                        stress,
                                                                                        lbmParams->GetStressParameter());
              }

              propertyCache.shearStressCache.Put(site.GetIndex(), stress);
            }

            if (propertyCache.vonMisesStressCache.RequiresRefresh())
            {
              distribn_t stress;
              StreamerImpl::CollisionType::CKernel::LatticeType::CalculateVonMisesStress(f_neq,
                                                                                         stress,
                                                                                         lbmParams->GetStressParameter());

              propertyCache.vonMisesStressCache.Put(site.GetIndex(), stress);
            }

            if (propertyCache.shearRateCache.RequiresRefresh())
            {
              distribn_t shear_rate = StreamerImpl::CollisionType::CKernel::LatticeType::CalculateShearRate(tau,
                                                                                                            f_neq,
                                                                                                            density);

              propertyCache.shearRateCache.Put(site.GetIndex(), shear_rate);
            }
          }

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BASESTREAMER_H */
