#ifndef HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
#define HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H

#include "mpiInclude.h"
#include "vis/DomainStats.h"
#include "vis/rayTracer/HSLToRGBConverter.h"
#include "vis/rayTracer/RayData.h"
#include "vis/VisSettings.h"
#include "vis/Vector3D.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      namespace DepthCuing
      {
	enum DepthCuing 
	{ 
	  NONE,
	  MIST,
	  DARKNESS
	};
      }

      template <DepthCuing::DepthCuing depthCuing>
	class RayDataEnhanced : public RayData<RayDataEnhanced<depthCuing> >
      {
      public:
      RayDataEnhanced() :
	mVelocitySum(0.0F),
	  mStressSum(0.0F),
	  mSurfaceNormalLightness(1.0F)
	  {
	  }
	
	//Processes the ray data for a normal (non wall) fluid site
	void DoUpdateDataForNormalFluidSite(const SiteData_t& iSiteData, 
					    const Vector3D<float>& iRayDirection,
					    const float iRayLengthInVoxel,
					    const DomainStats& iDomainStats,
					    const VisSettings& iVisSettings)
	{
	  mVelocitySum += 
	    iSiteData.Velocity * 
	    iRayLengthInVoxel *
	    (float) iDomainStats.velocity_threshold_max_inv;
	
	  if (iVisSettings.mStressType == lb::VonMises)
	  {
	    // update the volume rendering of the von Mises stress flow field
	    mStressSum = iSiteData.Stress *
	      iRayLengthInVoxel *
	      (float) iDomainStats.stress_threshold_max_inv;
	  }
	}
	
	void DoUpdateDataForWallSite(const SiteData_t& iSiteData, 
				     const Vector3D<float>& iRayDirection,
				     const float iRayLengthInVoxel,
				     const DomainStats& iDomainStats,
				     const VisSettings& iVisSettings,
				     const double* iWallNormal)
	{
	  DoUpdateDataForNormalFluidSite(iSiteData,
					 iRayDirection,
					 iRayLengthInVoxel,
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

	  //Scale the surface normal lightness between mParallelSurfaceAttenuation
	  //and 1.0F
	  mSurfaceNormalLightness *= (mParallelSurfaceAttenuation + 
				      (1.0F - mParallelSurfaceAttenuation)*fabs(lDotProduct)); 
	}
	  
	void DoGetVelocityColour(unsigned char oColour[3],
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
 
	void DoGetStressColour(unsigned char oColour[3],
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

	void DoMergeIn(const RayDataEnhanced& iOtherRayData,
		       const VisSettings& iVisSettings)
	{
	  mVelocitySum += iOtherRayData.GetVelocitySum();

	  if (iVisSettings.mStressType != lb::ShearStress)
	  {
	    mStressSum += iOtherRayData.GetStressSum();;
	  }
	
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
	  return GetVelocitySum()/
	    RayData<RayDataEnhanced<depthCuing> >::GetCumulativeLengthInFluid();
	}
      
	float GetAverageStress() const
	{
	  return GetStressSum()/
	    RayData<RayDataEnhanced<depthCuing> >::GetCumulativeLengthInFluid();
	}

      private:
	float GetLightnessValue(const float iNormalisedDistance) const
	{
	  assert(GetSurfaceNormalLightness() >= 0.0F && GetSurfaceNormalLightness() <= 1.0F);

	  if(depthCuing == DepthCuing::NONE)
	  {
	    return mLowestLightness +
	      GetSurfaceNormalLightness()*mSurfaceNormalLightnessRange;
	  }
	  else if (depthCuing == DepthCuing::MIST)
	  {
	    //Set the smallest lightness value to between
            //the mimimum lightness and 1.0F based on the normalised distance between
	    //the viewpoint and the first cluster hit
	    float lLightnessValue = mLowestLightness +
	      GetSurfaceNormalLightness()*mSurfaceNormalLightnessRange +
	      (1.0F - mLowestLightness)*iNormalisedDistance;

	    if(lLightnessValue > 1.0F)
	    {
	      return 1.0F;
	    } 
	    return lLightnessValue;
	  }
	  else if (depthCuing == DepthCuing::DARKNESS)
	  {
	    //Set the maximum lightness to be between 0.8F and mLowestLighness
	    //based on the noramlised distance
	    float lLightnessValue = (GetSurfaceNormalLightness() - 1.0F) *
	      mSurfaceNormalLightnessRange +
	      0.8F*(1.0F - iNormalisedDistance);
		    
	    if (lLightnessValue < mLowestLightness)
	    {
	      return mLowestLightness;
	    }
	    return lLightnessValue;
	  }
	  
	}

	float mSurfaceNormalLightness;
	float mVelocitySum;
	float mStressSum;

	
	const static float mSurfaceNormalLightnessRange;
	const static float mParallelSurfaceAttenuation;

	const static float mLowestLightness;
      };
    }
  }

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
