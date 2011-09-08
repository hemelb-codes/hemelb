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
      const float RayDataEnhanced::mIntensityMultipleThroughPerpendicularWalls = 0.75F;
      const float RayDataEnhanced::mLowestLightness = 0.3F;
      const float RayDataEnhanced::mHighestLightness = 0.8F;


      RayDataEnhanced::RayDataEnhanced() :
	mVelocitySum(0.0F),
	mStressSum(0.0F),
	mIntensity(1.0F)
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

	//We want the attentuation to be between
	//mIntensityMultipleThroughPerpendicularWalls and 1
	float lIntensityMultiple = mIntensityMultipleThroughPerpendicularWalls +
	  (1.0F-mIntensityMultipleThroughPerpendicularWalls)*lDotProduct;

	//Scale the current intensity of the ray
	mIntensity *= lIntensityMultiple; 
      }
      
      void RayDataEnhanced::DoGetVelocityColour(unsigned char oColour[3]) const
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
					   GetLightnessValue(),
					   oColour);
      }
    

      void RayDataEnhanced::DoGetStressColour(unsigned char oColour[3]) const
      {
	float lStressSaturation =
	  util::NumericalFunctions::enforceBounds<float>(0.5F+7.5F*GetAverageStress(),
							 0.0F,
							 1.0F);
	
	HSLToRGBConverter::ConvertHSLToRGB(230.0F,
					   lStressSaturation,
					   GetLightnessValue(),
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
	
	mIntensity *= iOtherRayData.GetIntensity();
      }
    
      float RayDataEnhanced::GetVelocitySum() const
      {
	return mVelocitySum;
      }
      
      float RayDataEnhanced::GetStressSum() const
      {
	return mStressSum;
      }
      
      float RayDataEnhanced::GetIntensity() const
      {
	return mIntensity;
      }

      float RayDataEnhanced::GetAverageVelocity() const
      {
	return GetVelocitySum()/GetCumulativeLengthInFluid();
      }
      
      float RayDataEnhanced::GetAverageStress() const
      {
	return GetStressSum()/GetCumulativeLengthInFluid();
      }

      float RayDataEnhanced::GetLightnessValue() const
      {
	assert(GetIntensity() >= 0.0F && mIntensity <= 1.0F);
	float lLightnessValue = GetIntensity()*
	  (mHighestLightness - mLowestLightness) +
	  mLowestLightness;
	
	assert(lLightnessValue >= mLowestLightness && mHighestLightness <= 1.0F);
	return lLightnessValue;
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
