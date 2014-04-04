// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
#define HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H

#include "net/mpi.h"
#include "vis/DomainStats.h"
#include "vis/rayTracer/HSLToRGBConverter.h"
#include "vis/rayTracer/RayData.h"
#include "vis/VisSettings.h"
#include "util/Vector3D.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      namespace DepthCuing
      {
        class Mist
        {
            static const float SurfaceNormalLightnessRange;
            static const float ParallelSurfaceAttenuation;
            static const float LowestLightness;
            static const float VelocityHueMin;
            static const float VelocityHueRange;
            static const float VelocitySaturation;
            static const float StressHue;
            static const float StressSaturationRange;
            static const float StressSaturationMin;

            //Obtain the lightness value for the ray, based on
            //the lightness obtained through surface normals
            //and optional depth cuing
            static float GetLightnessValue(const float normalisedDistance,
                                           const float surfaceNormalLightness)
            {
              //Set the smallest lightness value to between
              //the mimimum lightness and 1.0F based on the normalised distance between
              //the viewpoint and the first cluster hit
              //Add onto this the surface normal lightness
              float lightnessValue = LowestLightness + (1.0F - LowestLightness)
                  * normalisedDistance + surfaceNormalLightness * SurfaceNormalLightnessRange;

              if (lightnessValue > 1.0F)
              {
                return 1.0F;
              }
              return lightnessValue;
            }
        };

        class None
        {
            static const float SurfaceNormalLightnessRange;
            static const float ParallelSurfaceAttenuation;
            static const float LowestLightness;
            static const float VelocityHueMin;
            static const float VelocityHueRange;
            static const float VelocitySaturation;
            static const float StressHue;
            static const float StressSaturationRange;
            static const float StressSaturationMin;

            static float GetLightnessValue(const float normalisedDistance,
                                           const float surfaceNormalLightness)
            {
              return LowestLightness + surfaceNormalLightness * SurfaceNormalLightnessRange;
            }
        };

        class Darkness
        {
            static const float SurfaceNormalLightnessRange;
            static const float ParallelSurfaceAttenuation;
            static const float LowestLightness;
            static const float VelocityHueMin;
            static const float VelocityHueRange;
            static const float VelocitySaturation;
            static const float StressHue;
            static const float StressSaturationRange;
            static const float StressSaturationMin;

            static float GetLightnessValue(const float normalisedDistance,
                                           const float surfaceNormalLightness)
            {
              //Set the maximum lightness to be between 0.8F and mLowestLighness
              //based on the normalised distance and take off the surface normal
              //lightness
              float lightnessValue = 0.8F * (1.0F - normalisedDistance) + (surfaceNormalLightness
                  - 1.0F) * SurfaceNormalLightnessRange;

              if (lightnessValue < LowestLightness)
              {
                return LowestLightness;
              }
              return lightnessValue;
            }
        };
      }

      /**
       * RayDataEnhanced - sums velocity and stress data of the ray trace and provides
       * with surface normal highlighting and optional depth cuing to enhance the
       * 3D perception.
       * NB functions prefixed Do should only be called by the base class
       */
      class RayDataEnhanced : public RayData<RayDataEnhanced>
      {
        public:
          RayDataEnhanced(int i, int j) :
            RayData<RayDataEnhanced> (i, j), mSurfaceNormalLightness(1.0F), mVelocitySum(0.0F),
                mStressSum(0.0F)
          {
          }

          RayDataEnhanced()
          {

          }

          //Processes the ray data for a normal (non wall) fluid site
          void DoUpdateDataForNormalFluidSite(const SiteData_t& iSiteData,
                                              const util::Vector3D<float>& iRayDirection,
                                              const float iRayLengthInVoxel,
                                              const VisSettings& iVisSettings)
          {
            //Add the velocity multiplied by the ray length in each voxel,
            mVelocitySum += iSiteData.velocity * iRayLengthInVoxel;

            if (iVisSettings.mStressType == lb::VonMises)
            {
              //Update the volume rendering of the von Mises stress flow field
              mStressSum = iSiteData.stress * iRayLengthInVoxel;
            }
          }

          //Processes the data for wall sites
          template<typename depthCuing>
          void DoUpdateDataForWallSite(const SiteData_t& iSiteData,
                                       const util::Vector3D<float>& iRayDirection,
                                       const float iRayLengthInVoxel,
                                       const VisSettings& iVisSettings,
                                       const util::Vector3D<double>* iWallNormal)
          {
            //Do everything that would be done for a normal fluid site
            DoUpdateDataForNormalFluidSite(iSiteData,
                                           iRayDirection,
                                           iRayLengthInVoxel,
                                           iVisSettings);

            double lDotProduct = iRayDirection.Dot(*iWallNormal);

            // Scale the surface normal lightness between mParallelSurfaceAttenuation
            // and 1.0F
            // Keep a copy for the special case
            mLastSurfaceNormalLightnessMultiplier = (depthCuing::ParallelSurfaceAttenuation + (1.0F
                - depthCuing::ParallelSurfaceAttenuation) * fabs(lDotProduct));

            mSurfaceNormalLightness *= mLastSurfaceNormalLightnessMultiplier;
          }

          void DoProcessTangentingVessel()
          {
            mSurfaceNormalLightness *= mLastSurfaceNormalLightnessMultiplier;
          }

          //Obtains the colour representing the velocity ray trace

          template<typename depthCuing>
          void DoGetVelocityColour(unsigned char oColour[3],
                                   const float iNormalisedDistanceToFirstCluster,
                                   const DomainStats& iDomainStats) const
          {
            float lVelocityHue = depthCuing::VelocityHueRange * GetAverageVelocity()
                * (float) (iDomainStats.velocity_threshold_max_inv) + depthCuing::VelocityHueMin;

            if (lVelocityHue >= 360.0F)
            {
              lVelocityHue = 0.0F;
            }

            HSLToRGBConverter::Convert(lVelocityHue,
                                       depthCuing::VelocitySaturation,
                                       depthCuing::GetLightnessValue(iNormalisedDistanceToFirstCluster,
                                                                     GetSurfaceNormalLightness()),
                                       oColour);
          }

          //Obtains the colour representing the stress ray trace
          template<typename depthCuing>
          void DoGetStressColour(unsigned char oColour[3],
                                 const float iNormalisedDistanceToFirstCluster,
                                 const DomainStats& iDomainStats) const
          {
            float
                lStressSaturation =
                    util::NumericalFunctions::enforceBounds<float>(depthCuing::StressSaturationMin
                                                                       + depthCuing::StressSaturationRange
                                                                           * GetAverageStress()
                                                                           * (float) (iDomainStats.stress_threshold_max_inv),
                                                                   0.0F,
                                                                   1.0F);

            HSLToRGBConverter::Convert(depthCuing::StressHue,
                                       lStressSaturation,
                                       depthCuing::GetLightnessValue(iNormalisedDistanceToFirstCluster,
                                                                     GetSurfaceNormalLightness()),
                                       oColour);
          }

          //Carries out the merging of the ray data in this
          //inherited type, for different segments of the same ray
          void DoCombine(const RayDataEnhanced& iOtherRayData)
          {
            //Add together velocities and stress sums
            mVelocitySum += iOtherRayData.GetVelocitySum();
            mStressSum += iOtherRayData.GetStressSum();
            //Multiply the surface lightness
            mSurfaceNormalLightness *= iOtherRayData.GetSurfaceNormalLightness();
          }

          float GetVelocitySum() const
          {
            return mVelocitySum;
          }

          float GetStressSum() const
          {
            return mStressSum;
          }

          float GetSurfaceNormalLightness() const
          {
            return mSurfaceNormalLightness;
          }

          float GetAverageVelocity() const
          {
            return GetVelocitySum() / RayData<RayDataEnhanced>::GetCumulativeLengthInFluid();
          }

          float GetAverageStress() const
          {
            return GetStressSum() / RayData<RayDataEnhanced>::GetCumulativeLengthInFluid();
          }

          static MPI_Datatype GetMpiType()
          {
            HEMELB_MPI_TYPE_BEGIN(type, RayDataEnhanced, 9);

            HEMELB_MPI_TYPE_ADD_MEMBER(i);
            HEMELB_MPI_TYPE_ADD_MEMBER(j);
            HEMELB_MPI_TYPE_ADD_MEMBER(mLengthBeforeRayFirstCluster);
            HEMELB_MPI_TYPE_ADD_MEMBER(mCumulativeLengthInFluid);
            HEMELB_MPI_TYPE_ADD_MEMBER(mDensityAtNearestPoint);
            HEMELB_MPI_TYPE_ADD_MEMBER(mStressAtNearestPoint);
            HEMELB_MPI_TYPE_ADD_MEMBER(mSurfaceNormalLightness);
            HEMELB_MPI_TYPE_ADD_MEMBER(mVelocitySum);
            HEMELB_MPI_TYPE_ADD_MEMBER(mStressSum);

            HEMELB_MPI_TYPE_END(type, RayDataEnhanced);

            return type;
          }

        private:
          float mSurfaceNormalLightness;
          float mVelocitySum;
          float mStressSum;

          // The last surface normal lightness multiplier
          // is retained for the case of tangenting a vessel
          // ie only passing through walls sites
          // NB: Keep this last as it isn't sent over MPI
          float mLastSurfaceNormalLightnessMultiplier;
      };
    }
  }

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
