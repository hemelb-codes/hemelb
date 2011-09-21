#ifndef HEMELB_VIS_RAY_H
#define HEMELB_VIS_RAY_H

#include <cassert>
#include <iostream>
#include <limits>

#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "vis/rayTracer/SiteData.h"
#include "vis/rayTracer/RayDataNormal.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      template <typename RayDataType>
	class Ray
      {
      public:
      Ray(util::Vector3D<float> iDirection) :
	mInWall(false)
	{
	  iDirection.Normalise();
	  mDirection = iDirection;
	
	  mInverseDirection = util::Vector3D<float>
	    ( 1.0F/iDirection.x,
	      1.0F/iDirection.y,
	      1.0F/iDirection.z );

       	}
      
	util::Vector3D<float> GetDirection() const
	{
	  return mDirection;
	}

	util::Vector3D<float> GetInverseDirection() const
	{
	  return mInverseDirection;
	}

	bool XIncreasing() const
	{
	  return GetDirection().x > 0.0F;
	}

	bool YIncreasing() const
	{
	  return GetDirection().y > 0.0F;
	}

	bool ZIncreasing() const
	{
	  return GetDirection().z > 0.0F;
	}

	void UpdateDataForWallSite(const SiteData_t& iSiteData,
			const float iRayLengthInVoxel,
			const float iRayUnitsInCluster,
			const DomainStats& iDomainStats,
			const VisSettings& iVisSettings,
			const double* iWallNormal)
	{
	  if (!mInWall)
	  {
	    mRayData.UpdateDataForWallSite(iSiteData,
					 GetDirection(),
					 iRayLengthInVoxel,
					 iRayUnitsInCluster + mRayUnitsTraversedToCluster,
					 iDomainStats,
					 iVisSettings,
					 iWallNormal);
	    mInWall = true;
	  }
	  else
	  {
	    UpdateDataForNormalFluidSite(iSiteData,
					 iRayLengthInVoxel,
					 iRayUnitsInCluster,
					 iDomainStats,
					 iVisSettings);
	  }
	}
	
	void UpdateDataForNormalFluidSite(const SiteData_t& iSiteData,
					  const float iRayLengthInVoxel,
					  const float iRayUnitsInCluster,
					  const DomainStats& iDomainStats,
					  const VisSettings& iVisSettings)
	{
	  mInWall = false;
	  
	  mRayData.UpdateDataForNormalFluidSite(iSiteData,
						GetDirection(),
						iRayLengthInVoxel,
						iRayUnitsInCluster + mRayUnitsTraversedToCluster,
						iDomainStats,
						iVisSettings);
	}
	

	RayDataType GetRayData()
	{
	  return mRayData;
	}

	bool CollectedNoData()
	{
	  return !mRayData.ContainsRayData();
	}
	
	void SetRayLengthTraversedToCluster(float iRayUnitsTraversedToCluster)
	{
	  mRayUnitsTraversedToCluster = iRayUnitsTraversedToCluster;
	}
	

      private:
	util::Vector3D<float> mDirection;
	util::Vector3D<float> mInverseDirection;

	bool mInWall;

	float mRayUnitsTraversedToCluster;

	RayDataType mRayData;
      };	
    }
  }
}

#endif // HEMELB_VIS_RAY_H
