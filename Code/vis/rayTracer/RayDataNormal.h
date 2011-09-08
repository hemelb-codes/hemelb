#ifndef HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H
#define HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H

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
      class RayDataNormal : public RayData<RayDataNormal>
      {
      public:
	RayDataNormal();
	
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
	
	void DoMergeIn(const RayDataNormal& iOtherRayData,
		       const VisSettings& iVisSettings);

	void DoGetVelocityColour(unsigned char oColour[3]) const;
 
	void DoGetStressColour(unsigned char oColour[3]) const;

	bool DoContainsRayData() const;

	void UpdateVelocityColour(float iDt, const float iPalette[3]);

	void UpdateStressColour(float iDt, const float iPalette[3]);
      
      private:
	float mVelR, mVelG, mVelB;
	float mStressR, mStressG, mStressB;
      };	
    }
  }

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H
