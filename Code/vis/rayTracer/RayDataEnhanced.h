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
	
	void DoUpdateData(const SiteData_t& iSiteData, 
			  const double* iWallNormal,
			  const Vector3D<float>& iRayDirection,
			  const float iRayLengthInVoxel,
			  const float iAbsoluteDistanceFromViewpoint,
			  const DomainStats& iDomainStats,
			  const VisSettings& iVisSettings);
	//Passed as references since pointer can't be
	//meaningly transfered over MPI
	
	void DoMergeIn(const RayDataEnhanced& iOtherRayData,
		       const VisSettings& iVisSettings);

	void DoGetVelocityColour(unsigned char oColour[3]) const;
 
	void DoGetStressColour(unsigned char oColour[3]) const;

	bool DoContainsRayData() const;
	
	float GetVelocitySum() const;
     
	float GetStressSum() const;
      
	float GetIntensity() const;

	float GetAverageVelocity() const;
      
	float GetAverageStress() const;

      private:
	float mVelocitySum;
	float mStressSum;
	float mIntensity;
	
	const static float mMinLogIntensityMultiple;
        };
    }
  }

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATAENHANCED_H
