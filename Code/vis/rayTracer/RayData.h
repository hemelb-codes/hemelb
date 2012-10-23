// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_RAYTRACER_RAYDATA_H
#define HEMELB_VIS_RAYTRACER_RAYDATA_H

#include "constants.h"
#include "mpiInclude.h"
#include "vis/BasicPixel.h"
#include "vis/DomainStats.h"
#include "vis/VisSettings.h"
#include "vis/rayTracer/SiteData.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      /**
       * The abstract base type for RayData, using the CRTP pattern
       * Contains functionality common to all derived types
       */

      template<typename Derived>
      class RayData : public BasicPixel
      {
        public:
          RayData(int i, int j) :
            BasicPixel(i, j)
          {
            // A cheap way of indicating no ray data
            mCumulativeLengthInFluid = 0.0F;
          }

          RayData()
          {

          }

          // Used to process the ray data for a normal (non-wall) fluid site
          // Calls the method common to all ray data processing classes 
          // and then the derived class method
          void UpdateDataForNormalFluidSite(const raytracer::SiteData_t& iSiteData,
                                            const util::Vector3D<float>& iRayDirection,
                                            const float iRayLengthInVoxel,
                                            const float iAbsoluteDistanceFromViewpoint,
                                            const DomainStats& iDomainStats,
                                            const VisSettings& iVisSettings)
          {
            UpdateRayDataCommon(iSiteData,
                                iRayDirection,
                                iRayLengthInVoxel,
                                iAbsoluteDistanceFromViewpoint,
                                iDomainStats,
                                iVisSettings);

            ((Derived*) (this))->DoUpdateDataForNormalFluidSite(iSiteData,
                                                                iRayDirection,
                                                                iRayLengthInVoxel,
                                                                iVisSettings);
          }

          // Processes the ray data for a normal (non wall) fluid site
          void UpdateDataForWallSite(const raytracer::SiteData_t& iSiteData,
                                     const util::Vector3D<float>& iRayDirection,
                                     const float iRayLengthInVoxel,
                                     const float iAbsoluteDistanceFromViewpoint,
                                     const DomainStats& iDomainStats,
                                     const VisSettings& iVisSettings,
                                     const util::Vector3D<double>* iWallNormal)
          {
            UpdateRayDataCommon(iSiteData,
                                iRayDirection,
                                iRayLengthInVoxel,
                                iAbsoluteDistanceFromViewpoint,
                                iDomainStats,
                                iVisSettings);

            ((Derived*) (this))->DoUpdateDataForWallSite(iSiteData,
                                                         iRayDirection,
                                                         iRayLengthInVoxel,
                                                         iVisSettings,
                                                         iWallNormal);
          }

          void ProcessTangentingVessel()
          {
            ((Derived*) (this))->DoProcessTangentingVessel();
          }

          // Merges in the data from another segment of ray (from another core)
          void Combine(const Derived& iOtherRayData)
          {
            // Carry out the merging specific to the derived class
            ((Derived*) (this))->DoCombine(iOtherRayData);

            // Sum length in fluid
            SetCumulativeLengthInFluid(GetCumulativeLengthInFluid() + iOtherRayData.GetCumulativeLengthInFluid());

            // Update data relating to site nearest to viewpoint
            if (iOtherRayData.GetLengthBeforeRayFirstCluster() < this->GetLengthBeforeRayFirstCluster())
            {
              SetLengthBeforeRayFirstCluster(iOtherRayData.GetLengthBeforeRayFirstCluster());

              SetNearestDensity(iOtherRayData.GetNearestDensity());
              SetNearestStress(iOtherRayData.GetNearestStress());
            }
          }

          // Obtains the colour representing the velocity ray trace
          void GetVelocityColour(unsigned char oColour[3],
                                 const VisSettings& iVisSettings,
                                 const DomainStats& iDomainStats) const
          {
            ((const Derived*) (this))->DoGetVelocityColour(oColour,
                                                           GetLengthBeforeRayFirstCluster()
                                                               / iVisSettings.maximumDrawDistance,
                                                           iDomainStats);
          }

          // Obtains the colour representing the stress ray trace
          void GetStressColour(unsigned char oColour[3],
                               const VisSettings& iVisSettings,
                               const DomainStats& iDomainStats) const
          {
            ((const Derived*) (this))->DoGetStressColour(oColour,
                                                         GetLengthBeforeRayFirstCluster()
                                                             / iVisSettings.maximumDrawDistance,
                                                         iDomainStats);
          }

          // Whether or not the data contained in the instance has valid
          // ray data, or if it has just been constructed
          bool ContainsRayData() const
          {
            return (GetCumulativeLengthInFluid() != 0.0F);
          }

          float GetNearestStress() const
          {
            return mStressAtNearestPoint;
          }

          float GetNearestDensity() const
          {
            return mDensityAtNearestPoint;
          }

          float GetCumulativeLengthInFluid() const
          {
            return mCumulativeLengthInFluid;
          }

          float GetLengthBeforeRayFirstCluster() const
          {
            return mLengthBeforeRayFirstCluster;
          }

          void LogDebuggingInformation() const
          {
            log::Logger::Log<log::Trace, log::OnePerCore>("Ray data at (%i,%i) with "
                                                            "(lengthToFirstCluster, lengthInFluid, nearestDensity, nearest stress) = (%f, %f, %f, %f)",
                                                          GetI(),
                                                          GetJ(),
                                                          GetLengthBeforeRayFirstCluster(),
                                                          GetCumulativeLengthInFluid(),
                                                          GetNearestDensity(),
                                                          GetNearestStress());
          }

        protected:
          static const float mLongestDistanceInVoxelInverse;

          static void PickColour(float value, float colour[3])
          {
            colour[0] = util::NumericalFunctions::enforceBounds<float>(4.F * value - 2.F, 0.F, 1.F);
            colour[1]
                = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * (float) fabs(value - 0.5F), 0.F, 1.F);
            colour[2] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * value, 0.F, 1.F);
          }

          float mLengthBeforeRayFirstCluster;
          float mCumulativeLengthInFluid;

          float mDensityAtNearestPoint;
          float mStressAtNearestPoint;

        private:
          // Perform processing of the ray data common to all derived
          // types, ie cumulative length in the fluid, nearest stress
          // and density.
          void UpdateRayDataCommon(const raytracer::SiteData_t& iSiteData,
                                   const util::Vector3D<float>& iRayDirection,
                                   const float iRayLengthInVoxel,
                                   const float iAbsoluteDistanceFromViewpoint,
                                   const DomainStats& iDomainStats,
                                   const VisSettings& iVisSettings)
          {
            if (GetCumulativeLengthInFluid() == 0.0F)
            {
              SetLengthBeforeRayFirstCluster(iAbsoluteDistanceFromViewpoint);

              // Keep track of the density nearest to the viewpoint
              SetNearestDensity( (iSiteData.density - (float) iDomainStats.density_threshold_min)
                  * (float) iDomainStats.density_threshold_minmax_inv);

              if (iVisSettings.mStressType == lb::VonMises || iVisSettings.mStressType == lb::ShearStress)
              {
                // Keep track of the stress nearest to the viewpoint
                SetNearestStress(iSiteData.stress * (float) (iDomainStats.stress_threshold_max_inv));
              }
            }

            SetCumulativeLengthInFluid(GetCumulativeLengthInFluid() + iRayLengthInVoxel);
          }

          void SetNearestStress(float iStress)
          {
            mStressAtNearestPoint = iStress;
          }

          void SetNearestDensity(float iDensity)
          {
            mDensityAtNearestPoint = iDensity;
          }

          void SetCumulativeLengthInFluid(float iLength)
          {
            mCumulativeLengthInFluid = iLength;
          }

          void SetLengthBeforeRayFirstCluster(float iLength)
          {
            mLengthBeforeRayFirstCluster = iLength;
          }
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_RAYDATA_H
