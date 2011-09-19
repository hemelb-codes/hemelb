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
	  
	void DoGetVelocityColour(unsigned char oColour[3],
				 const float iNormalisedDistance) const;
 
	void DoGetStressColour(unsigned char oColour[3],
			       const float iNormalisedDistance) const;

	void DoMergeIn(const RayDataEnhanced& iOtherRayData,
		       const VisSettings& iVisSettings);

	bool DoContainsRayData() const;
	
	float GetVelocitySum() const;
     
	float GetStressSum() const;
      
	float GetSurfaceNormalLightness() const;

	float GetAverageVelocity() const;
      
	float GetAverageStress() const;

      private:
	float GetLightnessValue(const float iNormalisedDistance) const;

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
