// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <iostream>

#include "util/Vector3D.h"
#include "vis/DomainStats.h"
#include "vis/rayTracer/RayDataNormal.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {

      RayDataNormal::RayDataNormal(int i, int j) :
        RayData<RayDataNormal> (i, j)
      {
        mVelR = 0.0F;
        mVelG = 0.0F;
        mVelB = 0.0F;

        mStressR = 0.0F;
        mStressG = 0.0F;
        mStressB = 0.0F;
      }

      RayDataNormal::RayDataNormal()
      {

      }

      void RayDataNormal::DoUpdateDataForNormalFluidSite(
                                                         const SiteData_t& iSiteData,
                                                         const util::Vector3D<float>& iRayDirection,
                                                         const float iRayLengthInVoxel,
                                                         const VisSettings& iVisSettings)
      {
        float lPalette[3];

        // update the volume rendering of the velocity flow field
        PickColour(iSiteData.velocity * (float) mDomainStats->velocity_threshold_max_inv, lPalette);

        UpdateVelocityColour(iRayLengthInVoxel, lPalette);

        if (iVisSettings.mStressType != lb::ShearStress)
        {
          // update the volume rendering of the von Mises stress flow field
          float lScaledStress = iSiteData.stress * (float) mDomainStats->stress_threshold_max_inv;

          PickColour(lScaledStress, lPalette);

          UpdateStressColour(iRayLengthInVoxel, lPalette);
        }
      }

      void RayDataNormal::DoUpdateDataForWallSite(const SiteData_t& iSiteData,
                                                  const util::Vector3D<float>& iRayDirection,
                                                  const float iRayLengthInVoxel,
                                                  const VisSettings& iVisSettings,
                                                  const util::Vector3D<double>* iWallNormal)
      {
        DoUpdateDataForNormalFluidSite(iSiteData, iRayDirection, iRayLengthInVoxel, iVisSettings);
      }

      void RayDataNormal::DoGetVelocityColour(unsigned char oColour[3],
                                              const float iNormalisedDistanceToFirstCluster,
                                              const DomainStats& iDomainStats) const
      {
        MakeColourComponent(mVelR * 255.0F, oColour[0]);
        MakeColourComponent(mVelG * 255.0F, oColour[1]);
        MakeColourComponent(mVelB * 255.0F, oColour[2]);
      }

      void RayDataNormal::DoGetStressColour(unsigned char oColour[3],
                                            const float iNormalisedDistanceToFirstCluster,
                                            const DomainStats& iDomainStats) const
      {
        MakeColourComponent(mStressR, oColour[0]);
        MakeColourComponent(mStressG, oColour[1]);
        MakeColourComponent(mStressB, oColour[2]);
      }

      void RayDataNormal::MakeColourComponent(float value, unsigned char& colour) const
      {
        colour
            = util::NumericalFunctions::enforceBounds<unsigned char>((unsigned char) (value
                                                                         / GetCumulativeLengthInFluid()),
                                                                     0,
                                                                     255);
      }

      void RayDataNormal::DoCombine(const RayDataNormal& iOtherRayData)
      {
        mVelR += iOtherRayData.mVelR;
        mVelG += iOtherRayData.mVelG;
        mVelB += iOtherRayData.mVelB;

        mStressR += iOtherRayData.mStressR;
        mStressG += iOtherRayData.mStressG;
        mStressB += iOtherRayData.mStressB;
      }

      void RayDataNormal::UpdateVelocityColour(float iDt, const float iPalette[3])
      {
        mVelR += iDt * iPalette[0];
        mVelG += iDt * iPalette[1];
        mVelB += iDt * iPalette[2];
      }

      void RayDataNormal::UpdateStressColour(float iDt, const float iPalette[3])
      {
        mStressR += iDt * iPalette[0];
        mStressG += iDt * iPalette[1];
        mStressB += iDt * iPalette[2];
      }

      void RayDataNormal::DoProcessTangentingVessel()
      {
      }

      const DomainStats* RayDataNormal::mDomainStats = NULL;
    }
  }

  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataNormal>::RegisterMpiDataType()
    {
      MPI_Datatype ret = vis::raytracer::RayDataNormal::GetMPIType();
      HEMELB_MPI_CALL(MPI_Type_commit, (&ret));
      return ret;
    }
  }
}

