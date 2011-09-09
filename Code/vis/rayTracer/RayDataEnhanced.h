#ifndef HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
#define HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H

#include "mpiInclude.h"
#include "vis/DomainStats.h"
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
      class RayDataEnhanced : public RayData<RayDataEnhanced>
      {
      public:
	RayDataEnhanced();
	
	void DoUpdateDataForNormalFluidSite(const SiteData_t& iSiteData, 
					    const Vector3D<float>& iRayDirection,
					    const float iRayLengthInVoxel,
					    const float iAbsoluteDistanceFromViewpoint,
					    const DomainStats& iDomainStats,
					    const VisSettings& iVisSettings);
	//Passed as references since pointer can't be
	//meaningly transfered over MPI
	
	void DoUpdateDataForWallSite(const SiteData_t& iSiteData, 
				     const Vector3D<float>& iRayDirection,
				     const float iRayLengthInVoxel,
				     const float iAbsoluteDistanceFromViewpoint,
				     const DomainStats& iDomainStats,
				     const VisSettings& iVisSettings,
				     const double* iWallNormal);
	  
	void DoGetVelocityColour(unsigned char oColour[3]) const;
 
	void DoGetStressColour(unsigned char oColour[3]) const;

	void DoMergeIn(const RayDataEnhanced& iOtherRayData,
		       const VisSettings& iVisSettings);

	bool DoContainsRayData() const;
	
	float GetVelocitySum() const;
     
	float GetStressSum() const;
      
	float GetIntensity() const;

	float GetAverageVelocity() const;
      
	float GetAverageStress() const;

	bool DoIsRayCompletelyAttenuated() const;

      private:
	float GetLightnessValue() const;
	
	void PerformDepthCuing
	  ( float iAbsoluteDistanceFromViewpoint,
	    const VisSettings& iVisSettings );

	float mIntensity;
	float mVelocitySum;
	float mStressSum;

	const static float mIntensityMultipleThroughPerpendicularWalls;

	const static float mLowestLightness;
	const static float mHighestLightness;
        };
    }
  }

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
