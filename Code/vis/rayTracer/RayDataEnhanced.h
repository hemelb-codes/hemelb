#ifndef HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
#define HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H

#include "mpiInclude.h"
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
	enum DepthCuing 
	{ 
	  NONE,
	  MIST,
	  DARKNESS
	};
      }

      /**
       * RayDataEnhanced - sums velocity and stress data of the ray trace and provides
       * with surface normal highlighting and optional depth cuing to enhance the
       * 3D perception.
       * NB functions prefixed Do should only be called by the base class
       */
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
					    const util::Vector3D<float>& iRayDirection,
					    const float iRayLengthInVoxel,
					    const VisSettings& iVisSettings)
	{
	  //Add the velocity multiplied by the ray length in each voxel,
	  mVelocitySum += 
	    iSiteData.GetVelocity() * 
	    iRayLengthInVoxel;
	
	  if (iVisSettings.mStressType == lb::VonMises)
	  {
	    //Update the volume rendering of the von Mises stress flow field
	    mStressSum = iSiteData.GetStress() *
	      iRayLengthInVoxel;
	  }
	}
	
	//Processes the data for wall sites
	void DoUpdateDataForWallSite(const SiteData_t& iSiteData, 
				     const util::Vector3D<float>& iRayDirection,
				     const float iRayLengthInVoxel,
				     const VisSettings& iVisSettings,
				     const double* iWallNormal)
	{ 
	  //Do everything that would be done for a normal fluid site
	  DoUpdateDataForNormalFluidSite(iSiteData,
					 iRayDirection,
					 iRayLengthInVoxel,
					 iVisSettings);

	  //Calculate the absolute dot product of the wall
	  //vector normal and the ray direction
	  util::Vector3D<float> lWallNormal = 
	    util::Vector3D<float>(static_cast<float>(iWallNormal[0]),
				  static_cast<float>(iWallNormal[1]),
				  static_cast<float>(iWallNormal[2]));
	  
	  float lDotProduct =
	   iRayDirection.DotProduct(lWallNormal);

	  // Scale the surface normal lightness between mParallelSurfaceAttenuation
	  // and 1.0F
	  // Keep a copy for the special case
	    mLastSurfaceNormalLightnessMultiplier = (mParallelSurfaceAttenuation + 
						     (1.0F - mParallelSurfaceAttenuation)*fabs(lDotProduct)); 

	    mSurfaceNormalLightness *= mLastSurfaceNormalLightnessMultiplier;
	}

	void DoProcessTangentingVessel()
	{
	  mSurfaceNormalLightness *= mLastSurfaceNormalLightnessMultiplier;
	}
	 
       	//Obtains the colour representing the velocity ray trace
	void DoGetVelocityColour(unsigned char oColour[3],
				 const float iNormalisedDistanceToFirstCluster,
	                         const DomainStats& iDomainStats) const
	{
	  float lVelocityHue = mVelocityHueRange*
	    GetAverageVelocity()*
	    static_cast<float>(iDomainStats.velocity_threshold_max_inv) +
	    mVelocityHueMin;

	  if (lVelocityHue >= 360.0F)
	  {
	    lVelocityHue = 0.0F;
	  }
	
	  HSLToRGBConverter::Convert(lVelocityHue,
				     mVelocitySaturation,
				     GetLightnessValue(iNormalisedDistanceToFirstCluster),
				     oColour);
	}
 
	//Obtains the colour representing the stress ray trace
	void DoGetStressColour(unsigned char oColour[3],
			       const float iNormalisedDistanceToFirstCluster,
			       const DomainStats& iDomainStats) const
	{
	  float lStressSaturation =
	    util::NumericalFunctions::enforceBounds<float>
	    ( mStressSaturationMin +
	      mStressSaturationRange*
	      GetAverageStress()*
	      static_cast<float>(iDomainStats.stress_threshold_max_inv),
	      0.0F,
	      1.0F );
	
	  HSLToRGBConverter::Convert(mStressHue,
				     lStressSaturation,
				     GetLightnessValue(iNormalisedDistanceToFirstCluster),
				     oColour);
	}

        //Carries out the merging of the ray data in this
	//inherited type, for different segments of the same ray
	void DoMergeIn(const RayDataEnhanced& iOtherRayData,
		       const VisSettings& iVisSettings)
	{
	  //Add together velocities and stress sums
	  mVelocitySum += iOtherRayData.GetVelocitySum();

	  if (iVisSettings.mStressType != lb::ShearStress)
	  {
	    mStressSum += iOtherRayData.GetStressSum();;
	  }
	
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
	  return GetVelocitySum()/
	    RayData<RayDataEnhanced<depthCuing> >::GetCumulativeLengthInFluid();
	}
      
	float GetAverageStress() const
	{
	  return GetStressSum()/
	    RayData<RayDataEnhanced<depthCuing> >::GetCumulativeLengthInFluid();
	}

      private:
	//Obtain the lightness value for the ray, based on
	//the lightness obtained through surface normals
	//and optional depth cuing
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
	    //Add onto this the surface normal lightness
	    float lLightnessValue = mLowestLightness +
	      (1.0F - mLowestLightness)*iNormalisedDistance +
	      GetSurfaceNormalLightness()*mSurfaceNormalLightnessRange;

	    if(lLightnessValue > 1.0F)
	    {
	      return 1.0F;
	    } 
	    return lLightnessValue;
	  }
	  else if (depthCuing == DepthCuing::DARKNESS)
	  {
	    //Set the maximum lightness to be between 0.8F and mLowestLighness
	    //based on the normalised distance and take off the surface normal
            //lightness
	    float lLightnessValue = 0.8F*(1.0F - iNormalisedDistance) +
	      (GetSurfaceNormalLightness() - 1.0F) *
	      mSurfaceNormalLightnessRange;
		    
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

	// The last surface normal lightness multiplier
	// is retained for the case of tangenting a vessel
	// ie only passing through walls sites
	// NB: Keep this last as it isn't sent over MPI
	float mLastSurfaceNormalLightnessMultiplier;
	
	const static float mSurfaceNormalLightnessRange;
	const static float mParallelSurfaceAttenuation;

	const static float mLowestLightness;

	const static float mVelocityHueMin;
	const static float mVelocityHueRange;
	const static float mVelocitySaturation;

	const static float mStressHue;
	const static float mStressSaturationRange;
	const static float mStressSaturationMin;	
      };
    }
  }

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
