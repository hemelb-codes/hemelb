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
#include "lb/kernels/BaseKernel.h"
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

        protected:
          template<bool tDoRayTracing, class LatticeType>
          inline static void UpdateMinsAndMaxes(const geometry::Site& site,
                                                const kernels::HydroVarsBase<LatticeType>& hydroVars,
                                                const LbmParameters* lbmParams,
                                                lb::MacroscopicPropertyCache& propertyCache)
          {
            if (propertyCache.densityCache.RequiresRefresh())
            {
              propertyCache.densityCache.Put(site.GetIndex(), hydroVars.density);
            }

            if (propertyCache.velocityCache.RequiresRefresh())
            {
              propertyCache.velocityCache.Put(site.GetIndex(), hydroVars.momentum / hydroVars.density);
            }

            if (propertyCache.wallShearStressMagnitudeCache.RequiresRefresh())
            {
              distribn_t stress;

              if (!site.IsEdge())
              {
                stress = NO_VALUE;
              }
              else
              {
                LatticeType::CalculateWallShearStressMagnitude(hydroVars.density,
                                                               hydroVars.GetFNeq().f,
                                                               site.GetWallNormal(),
                                                               stress,
                                                               lbmParams->GetStressParameter());
              }

              propertyCache.wallShearStressMagnitudeCache.Put(site.GetIndex(), stress);
            }

            if (propertyCache.vonMisesStressCache.RequiresRefresh())
            {
              distribn_t stress;
              StreamerImpl::CollisionType::CKernel::LatticeType::CalculateVonMisesStress(hydroVars.GetFNeq().f,
                                                                                         stress,
                                                                                         lbmParams->GetStressParameter());

              propertyCache.vonMisesStressCache.Put(site.GetIndex(), stress);
            }

            if (propertyCache.shearRateCache.RequiresRefresh())
            {
              distribn_t shear_rate =
                  StreamerImpl::CollisionType::CKernel::LatticeType::CalculateShearRate(hydroVars.tau,
                                                                                        hydroVars.GetFNeq().f,
                                                                                        hydroVars.density);

              propertyCache.shearRateCache.Put(site.GetIndex(), shear_rate);
            }

            if (propertyCache.stressTensorCache.RequiresRefresh())
            {
              util::Matrix3D stressTensor;
              StreamerImpl::CollisionType::CKernel::LatticeType::CalculateStressTensor(hydroVars.density,
                                                                                       hydroVars.tau,
                                                                                       hydroVars.GetFNeq().f,
                                                                                       stressTensor);

              propertyCache.stressTensorCache.Put(site.GetIndex(), stressTensor);

            }

            if (propertyCache.tractionVectorCache.RequiresRefresh())
            {
              util::Vector3D<LatticeStress> tractionOnAPoint(0);

              /*
               * Wall normals are only available at the sites marked as being at the domain edge.
               * For the sites in the fluid bulk, the traction vector will be 0.
               */
              if (site.IsEdge())
              {
                LatticeType::CalculateTractionVectorOnAPoint(hydroVars.density,
                                                             hydroVars.tau,
                                                             hydroVars.GetFNeq().f,
                                                             site.GetWallNormal(),
                                                             tractionOnAPoint);
              }

              propertyCache.tractionVectorCache.Put(site.GetIndex(), tractionOnAPoint);

            }

            if (propertyCache.tangentialProjectionTractionVectorCache.RequiresRefresh())
            {
              util::Vector3D<LatticeStress> tangentialProjectionTractionOnAPoint(0);

              /*
               * Wall normals are only available at the sites marked as being at the domain edge.
               * For the sites in the fluid bulk, the traction vector will be 0.
               */
              if (site.IsEdge())
              {
                LatticeType::CalculateTangentialProjectionTractionVector(hydroVars.density,
                                                                         hydroVars.tau,
                                                                         hydroVars.GetFNeq().f,
                                                                         site.GetWallNormal(),
                                                                         tangentialProjectionTractionOnAPoint);
              }

              propertyCache.tangentialProjectionTractionVectorCache.Put(site.GetIndex(), tangentialProjectionTractionOnAPoint);

            }
          }
      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_BASESTREAMER_H */
