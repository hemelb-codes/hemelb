#ifndef HEMELB_VIS_RAYDATA_H
#define HEMELB_VIS_RAYDATA_H

//#define NDEBUG
#include "assert.h"


#include "constants.h"
#include "mpiInclude.h"
#include "vis/DomainStats.h"
#include "vis/VisSettings.h"
#include "vis/rayTracer/SiteData.h"
#include "lb/LbmParameters.h"
#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      template<typename Derived>
	class RayData
      {
      public:
	RayData()
	{
	  //A cheap way of indicating no ray data
	  mCumulativeLengthInFluid = 0.0F;
	
	  mStressAtNearestPoint = NO_VALUE_F;
	  mDensityAtNearestPoint = NO_VALUE_F;
	  mLengthBeforeRayFirstCluster = NO_VALUE_F;
	}

	void UpdateDataForNormalFluidSite(const raytracer::SiteData_t& iSiteData,
					  const Vector3D<float>& iRayDirection,
					  const float iRayLengthInVoxel,
					  const float iAbsoluteDistanceFromViewpoint,
					  const DomainStats& iDomainStats,
					  const VisSettings& iVisSettings)
	{
	  UpdateDataCommon(iSiteData,
			   iRayDirection,
			   iRayLengthInVoxel,
			   iAbsoluteDistanceFromViewpoint,
			   iDomainStats,
			   iVisSettings);

	  static_cast<Derived*>(this)->DoUpdateDataForNormalFluidSite
	    ( iSiteData, 
	      iRayDirection,
	      iRayLengthInVoxel,
	      iAbsoluteDistanceFromViewpoint,
	      iDomainStats,
	      iVisSettings );
	}

	void UpdateDataForWallSite(const raytracer::SiteData_t& iSiteData,
				   const Vector3D<float>& iRayDirection,
				   const float iRayLengthInVoxel,
				   const float iAbsoluteDistanceFromViewpoint,
				   const DomainStats& iDomainStats,
				   const VisSettings& iVisSettings,
				   const double* iWallNormal)
	{
	  assert(iWallNormal != NULL);

	  UpdateDataCommon(iSiteData,
			   iRayDirection,
			   iRayLengthInVoxel,
			   iAbsoluteDistanceFromViewpoint,
			   iDomainStats,
			   iVisSettings);

	  static_cast<Derived*>(this)->DoUpdateDataForWallSite
	    ( iSiteData, 
	      iRayDirection,
	      iRayLengthInVoxel,
	      iAbsoluteDistanceFromViewpoint,
	      iDomainStats,
	      iVisSettings,
	      iWallNormal);
	}

	void MergeIn(const Derived& iOtherRayData, const VisSettings& iVisSettings)
	{
	  static_cast<Derived*>(this)->DoMergeIn(iOtherRayData, iVisSettings);

	  ///Sum length in fluid
	  SetCumulativeLengthInFluid(
	    GetCumulativeLengthInFluid() + iOtherRayData.GetCumulativeLengthInFluid());

	  assert(GetCumulativeLengthInFluid() > 0.0F);

	  //Update data relating to site nearest to viewpoint
	  if (iOtherRayData.GetLengthBeforeRayFirstCluster() 
	      < this->GetLengthBeforeRayFirstCluster())
	  {
	    SetLengthBeforeRayFirstCluster(
	      iOtherRayData.GetLengthBeforeRayFirstCluster());
	    
	    SetNearestDensity(iOtherRayData.GetNearestDensity());
	    SetNearestStress(iOtherRayData.GetNearestStress());
	  }
	}

	void GetVelocityColour(unsigned char oColour[3]) const
	{
	  return static_cast<const Derived*>(this)->DoGetVelocityColour(oColour);
	}

	void GetStressColour(unsigned char oColour[3]) const
	{
	  return static_cast<const Derived*>(this)->DoGetStressColour(oColour);
	}

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

	bool IsRayCompletelyAttenuated() const
	{
	  return static_cast<const Derived*>(this)->DoIsRayCompletelyAttenuated();
	}

      private:
	float mLengthBeforeRayFirstCluster;
	float mCumulativeLengthInFluid;

	float mDensityAtNearestPoint;
	float mStressAtNearestPoint;
	
	void UpdateDataCommon(const raytracer::SiteData_t& iSiteData,
			      const Vector3D<float>& iRayDirection,
			      const float iRayLengthInVoxel,
			      const float iAbsoluteDistanceFromViewpoint,
			      const DomainStats& iDomainStats,
			      const VisSettings& iVisSettings)
	{
	  assert(iSiteData.Density != 0.0F);
	 
	  if (GetCumulativeLengthInFluid() == 0.0F)
	  {
	    SetLengthBeforeRayFirstCluster(iAbsoluteDistanceFromViewpoint);

	    // Keep track of the density nearest to the viewpoint
	    SetNearestDensity((iSiteData.Density - 
			       (float) iDomainStats.density_threshold_min) *
			      (float) iDomainStats.density_threshold_minmax_inv);
	  
	    if (iVisSettings.mStressType == lb::VonMises ||
		iVisSettings.mStressType == lb::ShearStress)
	    {
	      // Keep track of the stress nearest to the viewpoint
	      SetNearestStress(iSiteData.Stress * 
			       static_cast<float>(iDomainStats.stress_threshold_max_inv));		
	    }
	  }

	  SetCumulativeLengthInFluid(
	    GetCumulativeLengthInFluid() + iRayLengthInVoxel);
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

#endif // HEMELB_VIS_RAYTRACERDATA_H
