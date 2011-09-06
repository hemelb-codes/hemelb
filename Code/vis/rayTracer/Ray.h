#ifndef HEMELB_VIS_RAY_H
#define HEMELB_VIS_RAY_H

#include <limits>

#include "vis/rayTracer/SiteData.h"
#include "vis/rayTracer/RayDataNormal.h"
#include "vis/Vector3D.h"

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
	Ray(Vector3D<float> iDirection)
	{
	  iDirection.Normalise();
	  mDirection = iDirection;
	
	  mInverseDirection = Vector3D<float>
	    ( 1.0F/iDirection.x,
	      1.0F/iDirection.y,
	      1.0F/iDirection.z );
	}
      
	Vector3D<float> GetDirection() const
	{
	  return mDirection;
	}

	Vector3D<float> GetInverseDirection() const
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

	void UpdateData(const SiteData_t& iSiteData,
			const double* iWallNormal,
			const float iRayLengthInVoxel,
			const float iRayUnitsInCluster,
			const DomainStats& iDomainStats,
			const VisSettings& iVisSettings)
	{
	  mRayData.UpdateData(iSiteData,
			      iWallNormal,
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
	Vector3D<float> mDirection;
	Vector3D<float> mInverseDirection;

	float mRayUnitsTraversedToCluster;

	RayDataType mRayData;
      };	
    }
  }
}

#endif // HEMELB_VIS_RAY_H
