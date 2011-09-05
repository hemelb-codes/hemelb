#ifndef HEMELB_VIS_RAYDATA_H
#define HEMELB_VIS_RAYDATA_H

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
          void UpdateData(const raytracer::SiteData_t& iSiteData,
                          const double* iWallNormal,
                          const float iRayLengthInVoxel,
			  const float iAbsoluteDistanceFromViewpoint,
			  const DomainStats& iDomainStats,
			  const VisSettings& iVisSettings)
          {
            static_cast<Derived*>(this)->DoUpdateData(iSiteData,
						      iWallNormal, 
						      iRayLengthInVoxel,
						      iAbsoluteDistanceFromViewpoint,
						      iDomainStats,
						      iVisSettings);
          }

          void MergeIn(const Derived& iOtherRayData, const VisSettings& iVisSettings)
          {
            static_cast<Derived*>(this)->DoMergeIn(iOtherRayData, iVisSettings);
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
            return static_cast<const Derived*>(this)->DoContainsRayData();
          }

      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACERDATA_H
