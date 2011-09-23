#ifndef HEMELB_VIS_RAYTRACER_RAYDATA_H
#define HEMELB_VIS_RAYTRACER_RAYDATA_H

#include "assert.h"

#include "constants.h"
#include "mpiInclude.h"
#include "vis/DomainStats.h"
#include "vis/VisSettings.h"
#include "vis/rayTracer/SiteData.h"
#include "lb/LbmParameters.h"
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
      class RayData
      {
        public:
          RayData()
          {
            // A cheap way of indicating no ray data
            mCumulativeLengthInFluid = 0.0F;
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

            static_cast<Derived*>(this)->DoUpdateDataForNormalFluidSite(iSiteData,
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
                                     const double* iWallNormal)
          {
            assert(iWallNormal != NULL);

            UpdateRayDataCommon(iSiteData,
                                iRayDirection,
                                iRayLengthInVoxel,
                                iAbsoluteDistanceFromViewpoint,
                                iDomainStats,
                                iVisSettings);

            static_cast<Derived*>(this)->DoUpdateDataForWallSite(iSiteData,
                                                                 iRayDirection,
                                                                 iRayLengthInVoxel,
                                                                 iVisSettings,
                                                                 iWallNormal);
          }

          void ProcessTangentingVessel()
          {
            static_cast<Derived*>(this)->DoProcessTangentingVessel();
          }

          // Merges in the data from another segment of ray (from another core)
          void MergeIn(const Derived& iOtherRayData, const VisSettings& iVisSettings)
          {
            // Carry out the merging specific to the derived class
            static_cast<Derived*>(this)->DoMergeIn(iOtherRayData, iVisSettings);

            // Sum length in fluid
            SetCumulativeLengthInFluid(GetCumulativeLengthInFluid()
                + iOtherRayData.GetCumulativeLengthInFluid());

            // Update data relating to site nearest to viewpoint
            if (iOtherRayData.GetLengthBeforeRayFirstCluster()
                < this->GetLengthBeforeRayFirstCluster())
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
            static_cast<const Derived*>(this)->DoGetVelocityColour(oColour,
								   GetLengthBeforeRayFirstCluster()
								   / iVisSettings.maximumDrawDistance,
								   iDomainStats);
          }

          // Obtains the colour representing the stress ray trace
          void GetStressColour(unsigned char oColour[3],
                               const VisSettings& iVisSettings,
                               const DomainStats& iDomainStats) const
          {
            static_cast<const Derived*>(this)->DoGetStressColour(oColour,
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

        protected:
          static const float mLongestDistanceInVoxelInverse;

        private:
          float mLengthBeforeRayFirstCluster;
          float mCumulativeLengthInFluid;

          float mDensityAtNearestPoint;
          float mStressAtNearestPoint;

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
            assert(iSiteData.GetDensity() != 0.0F);

            if (GetCumulativeLengthInFluid() == 0.0F)
            {
              SetLengthBeforeRayFirstCluster(iAbsoluteDistanceFromViewpoint);

              // Keep track of the density nearest to the viewpoint
              SetNearestDensity( (iSiteData.GetDensity()
                  - (float) iDomainStats.density_threshold_min)
                  * (float) iDomainStats.density_threshold_minmax_inv);

              if (iVisSettings.mStressType == lb::VonMises
                  || iVisSettings.mStressType == lb::ShearStress)
              {
                // Keep track of the stress nearest to the viewpoint
                SetNearestStress(iSiteData.GetStress()
                    * static_cast<float>(iDomainStats.stress_threshold_max_inv));
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
