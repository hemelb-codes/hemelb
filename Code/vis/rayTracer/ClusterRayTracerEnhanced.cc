//#define NDEBUG;
#include <cassert>

#include "vis/rayTracer/ClusterRayTracerEnhanced.h"
#include "vis/rayTracer/RayEnhanced.h"

#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      ClusterRayTracerEnhanced::ClusterRayTracerEnhanced
      (	const Viewpoint& iViewpoint,
	Screen& iScreen,
	const DomainStats& iDomainStats,
	const VisSettings& iVisSettings,
	const hemelb::geometry::LatticeData& iLatticeData) :
	ClusterRayTracer(iViewpoint, 
			 iScreen, 
			 iDomainStats,
			 iVisSettings,
			 iLatticeData)
      {
      }

      void ClusterRayTracerEnhanced::CastRayForPixel
      (const Cluster& iCluster,
       const XYCoordinates<int>& iPixel,
       const Vector3D<float>& iRayDirection)
      {
	RayEnhanced lRay(iRayDirection);

	//These tell us how many ray units get us into the cluster
	//and after how many ray units we are out
	float lMaximumRayUnits;
	float lMinimumRayUnits;
	GetRayUnitsFromViewpointToCluster
	  (lRay, lMaximumRayUnits, lMinimumRayUnits);
	
	CastRay(iCluster, lRay, lMaximumRayUnits, lMinimumRayUnits);

	//Make sure the ray hasn't reached infinity
	if (lRay.LengthToFirstRayIntersection != std::numeric_limits<float>::max())
	{
	  ColPixel col_pixel(iPixel.x, iPixel.y, lRay.LengthToFirstRayIntersection + lMinimumRayUnits, lRay.Length, 
			     (lRay.Density - (float) mDomainStats.density_threshold_min)
			     * (float) mDomainStats.density_threshold_minmax_inv, lRay.Stress
			     != std::numeric_limits<float>::max()
			     ? lRay.Stress * (float) mDomainStats.stress_threshold_max_inv
			     : std::numeric_limits<float>::max(), lRay.VelocityColour, lRay.StressColour);

	  mScreen.AddPixel(&col_pixel, &mVisSettings);
	}
      }

      void ClusterRayTracer::UpdateRayData
      (const Cluster& iCluster,
       site_t iBlockNumber,
       site_t iSiteNumber,
       float iLengthFromClusterFirstIntersectionToVoxel,
       float iRayLengthInVoxel,
       Ray& ioRay)
      {
		const SiteData_t* lSiteData = iCluster.GetSiteData(iBlockNumber, iSiteNumber);
	
	if (lSiteData->Density < 0.0F)
	{
	  return; // solid voxel
	}

	float lPalette[3];

	// update the volume rendering of the velocity flow field
	ColPixel::PickColour(lSiteData->Velocity * (float) mDomainStats.velocity_threshold_max_inv,
			     lPalette);

	UpdateColour(iRayLengthInVoxel, lPalette, ioRay.VelocityColour);

	if (mVisSettings.mStressType != lb::ShearStress)
	{
	  // update the volume rendering of the von Mises stress flow field
	  float lScaledStress = lSiteData->Stress * (float) mDomainStats.stress_threshold_max_inv;

	  ColPixel::PickColour(lScaledStress, lPalette);

	  UpdateColour(iRayLengthInVoxel, lPalette, ioRay.StressColour);
	}

	ioRay.Length += iRayLengthInVoxel;

	if (ioRay.Density < 0.0F)
	{
	  ioRay.LengthToFirstRayIntersection = iLengthFromClusterFirstIntersectionToVoxel;

	  // keep track of the density nearest to the view point
	  ioRay.Density = lSiteData->Density;
	  
	  // keep track of the stress nearest to the view point
	  ioRay.Stress = lSiteData->Stress;	
	}
      }

    }
  }
}
