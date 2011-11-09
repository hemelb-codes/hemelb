#ifndef HEMELB_VIS_RAY_H
#define HEMELB_VIS_RAY_H

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
      template<typename RayDataType>
      class Ray
      {
        public:
          Ray(util::Vector3D<float> iDirection, int i, int j) :
            mInWall(false), mPassedThroughNormalFluidSite(false), mRayData(i, j)
          {
            iDirection.Normalise();
            mDirection = iDirection;

            mInverseDirection = util::Vector3D<float>(1.0F / iDirection.x,
                                                      1.0F / iDirection.y,
                                                      1.0F / iDirection.z);

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
            //Have we just entered the wall?
            if (!mInWall)
            {
              mRayData.UpdateDataForWallSite(iSiteData,
                                             GetDirection(),
                                             iRayLengthInVoxel,
                                             iRayUnitsInCluster + mRayUnitsTraversedToCluster,
                                             iDomainStats,
                                             iVisSettings,
                                             iWallNormal);
              //We're in the wall
              mInWall = true;
            }
            else
            {
              //We've already processed a wall site - process as a normal site
              mRayData.UpdateDataForNormalFluidSite(iSiteData,
                                                    GetDirection(),
                                                    iRayLengthInVoxel,
                                                    iRayUnitsInCluster
                                                        + mRayUnitsTraversedToCluster,
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
            //Set mInWall to false in case we've just left a wall
            mInWall = false;

            //We know we've passed through normal fluid
            mPassedThroughNormalFluidSite = true;

            mRayData.UpdateDataForNormalFluidSite(iSiteData,
                                                  GetDirection(),
                                                  iRayLengthInVoxel,
                                                  iRayUnitsInCluster + mRayUnitsTraversedToCluster,
                                                  iDomainStats,
                                                  iVisSettings);
          }

          void ProcessSolidSite()
          {
            //Special case - tangenting the vessel and never reaching
            //normal fluid sites
            if (mInWall && !mPassedThroughNormalFluidSite)
            {
              mRayData.ProcessTangentingVessel();
            }

            //We're out the wall
            mInWall = false;
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

          //mInWall indicates whether the ray is in a wall or not 
          //- if a wall site has just been processed
          bool mInWall;

          //mPassedThroughNormalFluid indicates whether the ray has passed
          //through normal fluid since entering a wall. If this remains false
          //and the ray hits a non-fluid site, the ray has just tangented
          //the surface of the vessel
          bool mPassedThroughNormalFluidSite;

          float mRayUnitsTraversedToCluster;

          RayDataType mRayData;
      };
    }
  }
}

#endif // HEMELB_VIS_RAY_H
