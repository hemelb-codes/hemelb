//#define NDEBUG;
#include <assert.h>
#include <cmath>
#include <iostream>

#include "util/utilityFunctions.h"
#include "vis/rayTracer/RayDataEnhanced.h"
#include "vis/rayTracer/HSLToRGBConverter.h"
#include "vis/ColPixel.h"
#include "vis/DomainStats.h"
#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {   
      const float RayDataEnhanced::mSurfaceNormalLightnessRange = 0.3F;
      const float RayDataEnhanced::mParallelSurfaceAttenuation = 0.5F;

      const float RayDataEnhanced::mLowestLightness = 0.3F;

      RayDataEnhanced::RayDataEnhanced() :
	mVelocitySum(0.0F),
	mStressSum(0.0F),
	mSurfaceNormalLightness(1.0F)
      {
      }

      void RayDataEnhanced::DoUpdateDataForNormalFluidSite
      ( const SiteData_t& iSiteData, 
	const Vector3D<float>& iRayDirection,
	const float iRayLengthInVoxel,
	const float iAbsoluteDistanceFromViewpoint,
	const DomainStats& iDomainStats,
	const VisSettings& iVisSettings )
      {
	mVelocitySum += iSiteData.Velocity * (float) iDomainStats.velocity_threshold_max_inv;
	
	iSiteData.Velocity * (float) iDomainStats.velocity_threshold_max_inv;

	if (iVisSettings.mStressType == lb::VonMises)
	{
	  // update the volume rendering of the von Mises stress flow field
	  mStressSum = iSiteData.Stress * (float) iDomainStats.stress_threshold_max_inv;
	}
      }

      void RayDataEnhanced::DoUpdateDataForWallSite(const SiteData_t& iSiteData, 
						    const Vector3D<float>& iRayDirection,
						    const float iRayLengthInVoxel,
						    const float iAbsoluteDistanceFromViewpoint,
						    const DomainStats& iDomainStats,
						    const VisSettings& iVisSettings,
						    const double* iWallNormal)
      {
	DoUpdateDataForNormalFluidSite(iSiteData,
				       iRayDirection,
				       iRayLengthInVoxel,
				       iAbsoluteDistanceFromViewpoint,
				       iDomainStats,
				       iVisSettings);

	//Calculate the absolute dot product of the wall
	//vector normal and the ray direction
	Vector3D<float> lWallNormal = 
	  Vector3D<float>(static_cast<float>(iWallNormal[0]),
			  static_cast<float>(iWallNormal[1]),
			  static_cast<float>(iWallNormal[2]));
	  
	float lDotProduct =
	  fabs(iRayDirection.DotProduct(lWallNormal));

	//Scale the current intensity of the ray by half the dot product 
	mSurfaceNormalLightness *= (mParallelSurfaceAttenuation + 
				    mParallelSurfaceAttenuation*fabs(lDotProduct)); 
      }
      
      void RayDataEnhanced::DoGetVelocityColour(unsigned char oColour[3],
						const float iNormalisedDistanceToFirstCluster) const
      {
	// We want the velocity hue to be between 240 degress
	// and 0 degrees

	float lVelocityHue = 120.0F*GetAverageVelocity() + 240.0F;

	if (lVelocityHue >= 360.0F)
	{
	  lVelocityHue = 0.0F;
	}
	
	HSLToRGBConverter::ConvertHSLToRGB(lVelocityHue,
					   1.0F,
					   GetLightnessValue(iNormalisedDistanceToFirstCluster),
					   oColour);
      }
    

      void RayDataEnhanced::DoGetStressColour(unsigned char oColour[3],
					      const float iNormalisedDistanceToFirstCluster) const
      {
	float lStressSaturation =
	  util::NumericalFunctions::enforceBounds<float>(0.5F+7.5F*GetAverageStress(),
							 0.0F,
							 1.0F);
	
	HSLToRGBConverter::ConvertHSLToRGB(230.0F,
					   lStressSaturation,
					   GetLightnessValue(iNormalisedDistanceToFirstCluster),
					   oColour);
      }

     
      void RayDataEnhanced::DoMergeIn(const RayDataEnhanced& iOtherRayData,
				      const VisSettings& iVisSettings)
      {
	mVelocitySum += iOtherRayData.GetVelocitySum();

	if (iVisSettings.mStressType != lb::ShearStress)
	{
	  mStressSum += iOtherRayData.GetStressSum();;
	}
	
	mSurfaceNormalLightness *= iOtherRayData.GetSurfaceNormalLightness();
      }
    
      float RayDataEnhanced::GetVelocitySum() const
      {
	return mVelocitySum;
      }
      
      float RayDataEnhanced::GetStressSum() const
      {
	return mStressSum;
      }
      
      float RayDataEnhanced::GetSurfaceNormalLightness() const
      {
	return mSurfaceNormalLightness;
      }

      float RayDataEnhanced::GetAverageVelocity() const
      {
	return GetVelocitySum()/GetCumulativeLengthInFluid();
      }
      
      float RayDataEnhanced::GetAverageStress() const
      {
	return GetStressSum()/GetCumulativeLengthInFluid();
      }

      float RayDataEnhanced::GetLightnessValue(const float iNormalisedDistance) const
      {
	assert(GetSurfaceNormalLightness() >= 0.0F && GetSurfaceNormalLightness() <= 1.0F);

	//To implement depth cuing, set the smallest lightness value to between 
        //the mimimum lightness and 1.0F based on the normalised distance between 
	//the viewpoint and the first cluster hit 
	float lLightnessValue = mLowestLightness + (1.0F - 0.3F)*iNormalisedDistance;
	
	//Map the surface normal lightness to between this value and
	//mSurfaceNormalLightnessRange above this
	lLightnessValue += GetSurfaceNormalLightness()*mSurfaceNormalLightnessRange;
	
	if (lLightnessValue < 1.0F)
	{
	  return lLightnessValue;
	}
	else 
	{
	  return 1.0F;
	}
      }
    }
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced>::RegisterMpiDataType()
  {
    int lRayDataEnhancedCount = 8;
    int lRayDataEnhancedBlocklengths[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype lRayDataEnhancedTypes[8] = { MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_UB };

    MPI_Aint lRayDataEnhancedDisps[8];

    lRayDataEnhancedDisps[0] = 0;

    for (int i = 1; i < lRayDataEnhancedCount; i++)
    {
      if (lRayDataEnhancedTypes[i - 1] == MPI_FLOAT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1]
	  + (sizeof(float) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_INT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(int) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_UNSIGNED)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(unsigned) * lRayDataEnhancedBlocklengths[i - 1]);
      }
    }
    MPI_Datatype type;
    MPI_Type_struct(lRayDataEnhancedCount,
		    lRayDataEnhancedBlocklengths,
		    lRayDataEnhancedDisps,
		    lRayDataEnhancedTypes,
		    &type);
    MPI_Type_commit(&type);
    return type;
  }
}
