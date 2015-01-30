// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_RAYTRACER_CLUSTERRAYTRACER_H
#define HEMELB_VIS_RAYTRACER_CLUSTERRAYTRACER_H

#include <cmath> 
#include <iostream>
#include <limits>

#include "geometry/SiteTraverser.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "vis/DomainStats.h"
#include "vis/PixelSet.h"
#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/ClusterTraverser.h"
#include "vis/rayTracer/Ray.h"
#include "vis/Screen.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      template<typename ClusterType, typename RayDataType>
      class ClusterRayTracer
      {
        public:
          ClusterRayTracer(const Viewpoint& iViewpoint,
                           Screen& iScreen,
                           const DomainStats& iDomainStats,
                           const VisSettings& iVisSettings,
                           const hemelb::geometry::LatticeData& iLatticeData,
                           const lb::MacroscopicPropertyCache& propertyCache) :
              viewpoint(iViewpoint), screen(iScreen), domainStats(iDomainStats), visSettings(iVisSettings), latticeData(iLatticeData), propertyCache(propertyCache)
          {
            // TODO: This is absolutely horrible, but neccessary until RayDataNormal is
            // removed. 
            RayDataNormal::mDomainStats = &iDomainStats;
          }

          void RenderCluster(const ClusterType& iCluster, PixelSet<RayDataType>& pixels)
          {
            mLowerSiteCordinatesOfClusterRelativeToViewpoint = iCluster.GetLeastSiteOnLeastBlockInImage()
                - viewpoint.GetViewpointLocation();

            //Calculate the projection of the cluster on the screen
            //refered to as the subimage
            CalculateSubImage(iCluster);

            //If the entire sub-image is off the screen,
            //no rendering is needed
            if (SubImageOffScreen())
            {
              return;
            }

            CropSubImageToScreen();

            CalculateVectorsToClusterSpanAndLowerLeftPixel(iCluster);

            CastRaysForEachPixel(iCluster, pixels);
          }

        private:
          void GetRayUnitsFromViewpointToCluster(const Ray<RayDataType> & iRay,
                                                 float & oMaximumRayUnits,
                                                 float & oMinimumRayUnits)
          {
            // (Remember that iRay.mDirection is normalised)
            float lMaxUnitRaysBasedOnX;
            float lMinUnitRaysBasedOnX;
            if (iRay.GetDirection().x > 0.0F)
            {
              lMaxUnitRaysBasedOnX = mViewpointCentreToMaxSite.x * iRay.GetInverseDirection().x;

              lMinUnitRaysBasedOnX = mViewpointCentreToMinSite.x * iRay.GetInverseDirection().x;
            }
            else if (iRay.GetDirection().x < 0.0F)
            {
              lMaxUnitRaysBasedOnX = mViewpointCentreToMinSite.x * iRay.GetInverseDirection().x;
              lMinUnitRaysBasedOnX = mViewpointCentreToMaxSite.x * iRay.GetInverseDirection().x;
            }
            else
            {
              lMaxUnitRaysBasedOnX = std::numeric_limits<float>::max();
              lMinUnitRaysBasedOnX = 0.0F;
            }

            float lMaxUnitRaysBasedOnY;
            float lMinUnitRaysBasedOnY;
            if (iRay.GetDirection().y > 0.0F)
            {
              lMaxUnitRaysBasedOnY = mViewpointCentreToMaxSite.y * iRay.GetInverseDirection().y;

              lMinUnitRaysBasedOnY = mViewpointCentreToMinSite.y * iRay.GetInverseDirection().y;
            }
            else if (iRay.GetDirection().y < 0.0F)
            {
              lMaxUnitRaysBasedOnY = mViewpointCentreToMinSite.y * iRay.GetInverseDirection().y;
              lMinUnitRaysBasedOnY = mViewpointCentreToMaxSite.y * iRay.GetInverseDirection().y;
            }
            else
            {
              lMaxUnitRaysBasedOnY = std::numeric_limits<float>::max();
              lMinUnitRaysBasedOnY = 0.0F;
            }

            float lMaxUnitRaysBasedOnZ;
            float lMinUnitRaysBasedOnZ;
            if (iRay.GetDirection().z > 0.0F)
            {
              lMaxUnitRaysBasedOnZ = mViewpointCentreToMaxSite.z * iRay.GetInverseDirection().z;

              lMinUnitRaysBasedOnZ = mViewpointCentreToMinSite.z * iRay.GetInverseDirection().z;
            }
            else if (iRay.GetDirection().z < 0.0F)
            {
              lMaxUnitRaysBasedOnZ = mViewpointCentreToMinSite.z * iRay.GetInverseDirection().z;
              lMinUnitRaysBasedOnZ = mViewpointCentreToMaxSite.z * iRay.GetInverseDirection().z;
            }
            else
            {
              lMaxUnitRaysBasedOnZ = std::numeric_limits<float>::max();
              lMinUnitRaysBasedOnZ = 0.0F;
            }

            //Maximum ray units from viewpoint to cluster
            //We want the minimum number - since at this point the ray is
            //completely out
            oMaximumRayUnits = std::min(std::min(lMaxUnitRaysBasedOnX, lMaxUnitRaysBasedOnY), lMaxUnitRaysBasedOnZ);

            //Maximum ray units to get us into the cluster
            //We want the maximum number - since only at this point
            // is the ray completely in
            oMinimumRayUnits = std::max(std::max(lMinUnitRaysBasedOnX, lMinUnitRaysBasedOnY), lMinUnitRaysBasedOnZ);

          }

          void CastRay(const ClusterType& iCluster,
                       Ray<RayDataType>& ioRay,
                       float iMaximumRayUnits,
                       float iMinimumRayUnits)
          {
            //It's possible for the ray to totally miss the cluster
            //This is because the sub-image is square while the cluster
            // projection won't be in most circumstances
            if (iMaximumRayUnits < iMinimumRayUnits)
            {
              return;
            }

            ioRay.SetRayLengthTraversedToCluster(iMinimumRayUnits);

            util::Vector3D<float> fromLowerSiteToFirstRayClusterIntersection = ioRay.GetDirection() * iMinimumRayUnits
                - mLowerSiteCordinatesOfClusterRelativeToViewpoint;

            TraverseBlocks(iCluster, fromLowerSiteToFirstRayClusterIntersection, ioRay);
          }

          void CalculateSubImage(const ClusterType& iCluster)
          {
            //The extent of the cluster when projected (ie the subimage)
            //is determined by projecting all eight vertices of the cuboid
            XYCoordinates<float> lSubImageLowerLeft = XYCoordinates<float>::MaxLimit();
            XYCoordinates<float> lSubImageUpperRight = XYCoordinates<float>::MinLimit();

            const std::vector<util::Vector3D<float> > lCorners = iCluster.GetCorners();

            for (std::vector<util::Vector3D<float> >::const_iterator lIt = lCorners.begin(); lIt != lCorners.end();
                lIt++)
            {
              UpdateSubImageExtentForCorner(*lIt, lSubImageLowerLeft, lSubImageUpperRight);
            }

            lowerLeftPixelCoordinatesOfSubImage =
                screen.template TransformScreenToPixelCoordinates<int>(lSubImageLowerLeft);

            // We add a unit vector here because the transformation will round down from float
            // to int.
            upperRightPixelCoordinatesOfSubImage =
                screen.template TransformScreenToPixelCoordinates<int>(lSubImageUpperRight);
          }

          void UpdateSubImageExtentForCorner(const util::Vector3D<float>& iCorner,
                                             XYCoordinates<float>& ioSubImageLowerLeft,
                                             XYCoordinates<float>& ioSubImageUpperRight)
          {
            XYCoordinates<float> lCornerProjection = viewpoint.FlatProject(iCorner);

            ioSubImageLowerLeft.UpdatePointwiseMin(lCornerProjection);
            ioSubImageUpperRight.UpdatePointwiseMax(lCornerProjection);
          }

          bool SubImageOffScreen()
          {
            return (lowerLeftPixelCoordinatesOfSubImage.x >= screen.GetPixelsX()
                || upperRightPixelCoordinatesOfSubImage.x < 0
                || lowerLeftPixelCoordinatesOfSubImage.y >= screen.GetPixelsY()
                || upperRightPixelCoordinatesOfSubImage.y < 0);
          }

          void CropSubImageToScreen()
          {
            lowerLeftPixelCoordinatesOfSubImage.x = util::NumericalFunctions::max(lowerLeftPixelCoordinatesOfSubImage.x,
                                                                                  0);

            upperRightPixelCoordinatesOfSubImage.x =
                util::NumericalFunctions::min(upperRightPixelCoordinatesOfSubImage.x, screen.GetPixelsX() - 1);

            lowerLeftPixelCoordinatesOfSubImage.y = util::NumericalFunctions::max(lowerLeftPixelCoordinatesOfSubImage.y,
                                                                                  0);

            upperRightPixelCoordinatesOfSubImage.y =
                util::NumericalFunctions::min(upperRightPixelCoordinatesOfSubImage.y, screen.GetPixelsY() - 1);
          }

          void CalculateVectorsToClusterSpanAndLowerLeftPixel(const ClusterType& iCluster)
          {
            mViewpointCentreToMaxSite = iCluster.GetMaxSite() - viewpoint.GetViewpointLocation();

            mViewpointCentreToMinSite = iCluster.GetMinSite() - viewpoint.GetViewpointLocation();

            //Obtaining the vector from the camera to the lower left pixel
            fromCameraToBottomLeftPixelOfSubImage = screen.GetCameraToBottomLeftOfScreenVector()
                + screen.GetPixelUnitVectorProjectionX() * (float) lowerLeftPixelCoordinatesOfSubImage.x
                + screen.GetPixelUnitVectorProjectionY() * (float) lowerLeftPixelCoordinatesOfSubImage.y;
          }

          void CastRaysForEachPixel(const ClusterType& iCluster, PixelSet<RayDataType>& pixels)
          {
            XYCoordinates<int> lPixel;

            //Loop over all the pixels
            util::Vector3D<float> lCameraToBottomRow = fromCameraToBottomLeftPixelOfSubImage;
            for (lPixel.x = lowerLeftPixelCoordinatesOfSubImage.x; lPixel.x <= upperRightPixelCoordinatesOfSubImage.x;
                ++lPixel.x)
            {
              util::Vector3D<float> lCameraToPixel = lCameraToBottomRow;
              for (lPixel.y = lowerLeftPixelCoordinatesOfSubImage.y; lPixel.y <= upperRightPixelCoordinatesOfSubImage.y;
                  ++lPixel.y)
              {
                CastRayForPixel(iCluster, lPixel, lCameraToPixel, pixels);

                lCameraToPixel += screen.GetPixelUnitVectorProjectionY();
              }

              lCameraToBottomRow += screen.GetPixelUnitVectorProjectionX();
            }
          }

          virtual void CastRayForPixel(const ClusterType& iCluster,
                                       const XYCoordinates<int>& iPixelCoordinates,
                                       const util::Vector3D<float>& iRayDirection,
                                       PixelSet<RayDataType>& pixels)
          {
            Ray<RayDataType> lRay(iRayDirection, iPixelCoordinates.x, iPixelCoordinates.y);

            //These tell us how many ray units get us into the cluster
            //and after how many ray units we are out
            float lMaximumRayUnits;
            float lMinimumRayUnits;
            GetRayUnitsFromViewpointToCluster(lRay, lMaximumRayUnits, lMinimumRayUnits);

            CastRay(iCluster, lRay, lMaximumRayUnits, lMinimumRayUnits);

            //Make sure the ray hasn't reached infinity
            if (!lRay.CollectedNoData())
            {
              pixels.AddPixel(lRay.GetRayData());
            }
          }

          void TraverseRayThroughBlock(const util::Vector3D<float>& fromFirstRayClusterIntersectionToLowerSiteOfCurrentBlock,
                                       const util::Vector3D<float>& iLocationInBlock,
                                       const ClusterType& iCluster,
                                       const util::Vector3D<site_t>& blockLocation,
                                       const site_t blockNumberOnCluster,
                                       float euclideanClusterLengthTraversedByRay,
                                       Ray<RayDataType>& ioRay)
          {
            //Work out which site we're currently in
            const util::Vector3D<site_t> truncatedLocationInBlock = RoundToNearestVoxel(iLocationInBlock);

            geometry::SiteTraverser siteTraverser(latticeData);
            siteTraverser.SetCurrentLocation(truncatedLocationInBlock);

            // In order to trace the rays through the voxels, we need to keep track of how far the
            // ray can travel to the next voxel in each of the three directions in ray units
            util::Vector3D<float> rayUnitsUntilNextSite =
                CalculateRayUnitsBeforeNextSite(fromFirstRayClusterIntersectionToLowerSiteOfCurrentBlock,
                                                util::Vector3D<float>(truncatedLocationInBlock),
                                                ioRay);

            while (siteTraverser.CurrentLocationValid())
            {
              // Firstly, work out in which direction we
              // can travel the least ray units before reaching
              // a vortex side
              const util::Direction::Direction directionOfLeastTravel = DirectionOfLeastTravel(rayUnitsUntilNextSite);

              // Find out how far the ray can move
              const float manhattanRayLengthThroughVoxel = rayUnitsUntilNextSite.GetByDirection(directionOfLeastTravel);

              const geometry::Block& block = latticeData.GetBlock(latticeData.GetBlockIdFromBlockCoords(blockLocation));

              if (!block.IsEmpty()) // Ensure fluid site
              {
                if (!block.SiteIsSolid(siteTraverser.GetCurrentIndex()))
                {
                  const site_t localContiguousId =
                      block.GetLocalContiguousIndexForSite(siteTraverser.GetCurrentIndex());

                  SiteData_t siteData;
                  siteData.density = propertyCache.densityCache.Get(localContiguousId);
                  siteData.velocity = propertyCache.velocityCache.Get(localContiguousId).GetMagnitude();

                  if (visSettings.mStressType == lb::ShearStress)
                  {
                    siteData.stress = propertyCache.wallShearStressMagnitudeCache.Get(localContiguousId);
                  }
                  else
                  {
                    siteData.stress = propertyCache.vonMisesStressCache.Get(localContiguousId);
                  }

                  const util::Vector3D<double>* lWallData = iCluster.GetWallData(blockNumberOnCluster,
                                                                                 siteTraverser.GetCurrentIndex());

                  if (lWallData == nullptr || lWallData->x == NO_VALUE)
                  {
                    ioRay.UpdateDataForNormalFluidSite(siteData,
                                                       manhattanRayLengthThroughVoxel
                                                           - euclideanClusterLengthTraversedByRay, // Manhattan Ray-length through the voxel
                                                       euclideanClusterLengthTraversedByRay, // euclidean ray units spent in cluster
                                                       domainStats,
                                                       visSettings);
                  }
                  else
                  {
                    ioRay.UpdateDataForWallSite(siteData,
                                                manhattanRayLengthThroughVoxel - euclideanClusterLengthTraversedByRay,
                                                euclideanClusterLengthTraversedByRay,
                                                domainStats,
                                                visSettings,
                                                lWallData);
                  }
                }
                else
                {
                  ioRay.ProcessSolidSite();
                }
              }
              else
              {
                ioRay.ProcessSolidSite();
              }

              //Update ray length traversed so far
              euclideanClusterLengthTraversedByRay = manhattanRayLengthThroughVoxel;

              //Update the block location and RayUnitsBeforeNextVoxel
              //in each direction
              switch (directionOfLeastTravel)
              {
                case util::Direction::X:
                  if (ioRay.XIncreasing())
                  {
                    siteTraverser.IncrementX();
                    rayUnitsUntilNextSite.x += ioRay.GetInverseDirection().x;
                  }
                  else
                  {
                    siteTraverser.DecrementX();
                    rayUnitsUntilNextSite.x -= ioRay.GetInverseDirection().x;
                  }

                  break;

                case util::Direction::Y:
                  if (ioRay.YIncreasing())
                  {
                    siteTraverser.IncrementY();
                    rayUnitsUntilNextSite.y += ioRay.GetInverseDirection().y;
                  }
                  else
                  {
                    siteTraverser.DecrementY();
                    rayUnitsUntilNextSite.y -= ioRay.GetInverseDirection().y;
                  }
                  break;

                case util::Direction::Z:
                  if (ioRay.ZIncreasing())
                  {
                    siteTraverser.IncrementZ();
                    rayUnitsUntilNextSite.z += ioRay.GetInverseDirection().z;
                  }
                  else
                  {
                    siteTraverser.DecrementZ();
                    rayUnitsUntilNextSite.z -= ioRay.GetInverseDirection().z;
                  }
                  break;
              }

            }
          }

          util::Vector3D<float> CalculateRayUnitsBeforeNextSite(const util::Vector3D<float>& iFirstRayClusterIntersectionToBlockLowerSite,
                                                                const util::Vector3D<site_t>& iTruncatedLocationInBlock,
                                                                const Ray<RayDataType>& iRay) const
          {
            util::Vector3D<float> lRayUnits;

            //The ray has already travelled iFirstRayClusterIntersectionToBlockLowerSite
            //If the ray is increasing in the co-ordinte it can travel as far
            //as the truncated location + 1, otherwise just the truncated location
            //for each co-ordinate

            //If the ray has zero in any direction, set the ray units to max

            if (iRay.GetDirection().x == 0.0F)
            {
              lRayUnits.x = std::numeric_limits<float>::max();
            }
            else
            {
              lRayUnits.x = iFirstRayClusterIntersectionToBlockLowerSite.x + (float) (iTruncatedLocationInBlock.x);
              if (iRay.XIncreasing())
              {
                lRayUnits.x += 1.0F;
              }
              //Convert from site units into ray units
              lRayUnits.x *= iRay.GetInverseDirection().x;
            }

            if (iRay.GetDirection().y == 0.0F)
            {
              lRayUnits.y = std::numeric_limits<float>::max();
            }
            else
            {
              lRayUnits.y = iFirstRayClusterIntersectionToBlockLowerSite.y + (float) (iTruncatedLocationInBlock.y);
              if (iRay.YIncreasing())
              {
                lRayUnits.y += 1.0F;
              }
              lRayUnits.y *= iRay.GetInverseDirection().y;
            }

            if (iRay.GetDirection().z == 0.0F)
            {
              lRayUnits.z = std::numeric_limits<float>::max();
            }
            else
            {
              lRayUnits.z = iFirstRayClusterIntersectionToBlockLowerSite.z + (float) (iTruncatedLocationInBlock.z);
              if (iRay.ZIncreasing())
              {
                lRayUnits.z += 1.0F;
              }
              lRayUnits.z *= iRay.GetInverseDirection().z;
            }

            return lRayUnits;
          }

          util::Vector3D<site_t> RoundToNearestVoxel(const util::Vector3D<float>& iUnboundLocation) const
          {
            util::Vector3D<site_t> lVoxelLocationInBlock;

            //Due to rounding errors, it's possible for the site location within a block
            //to be outside the wrong block
            lVoxelLocationInBlock.x = util::NumericalFunctions::enforceBounds<site_t>((site_t) iUnboundLocation.x,
                                                                                      0,
                                                                                      latticeData.GetBlockSize() - 1);

            lVoxelLocationInBlock.y = util::NumericalFunctions::enforceBounds<site_t>((site_t) iUnboundLocation.y,
                                                                                      0,
                                                                                      latticeData.GetBlockSize() - 1);

            lVoxelLocationInBlock.z = util::NumericalFunctions::enforceBounds<site_t>((site_t) iUnboundLocation.z,
                                                                                      0,
                                                                                      latticeData.GetBlockSize() - 1);

            return lVoxelLocationInBlock;
          }

          util::Direction::Direction DirectionOfLeastTravel(const util::Vector3D<float>& iRayUnitsBeforeNextVoxelOrBlock) const
          {

            if (iRayUnitsBeforeNextVoxelOrBlock.x < iRayUnitsBeforeNextVoxelOrBlock.y)
            {
              //X is less than Y
              if (iRayUnitsBeforeNextVoxelOrBlock.x < iRayUnitsBeforeNextVoxelOrBlock.z)
              {
                //X is less than Y and Z
                return util::Direction::X;
              }
              else
              {
                //X is less than Y
                //Z is less Than X (And Y)
                return util::Direction::Z;
              }
            }
            else
            {
              // Y is less than X
              if (iRayUnitsBeforeNextVoxelOrBlock.y < iRayUnitsBeforeNextVoxelOrBlock.z)
              {
                //Y is less than X and Z
                return util::Direction::Y;
              }
              else
              {
                //Y is less than X
                //Z is less than Y (and so X)
                return util::Direction::Z;
              }

            }
          }

          void TraverseBlocks(const ClusterType& iCluster,
                              const util::Vector3D<float>& fromLowerClusterSiteToFirstRayIntersection,
                              Ray<RayDataType>& ioRay)
          {
            float blockSizeAsFloat = (float) (latticeData.GetBlockSize());

            //Calculate the coordinates of the block within the cluster where
            //the ray first intersects
            const util::Vector3D<site_t> blockHoldingFirstIntersection =
                GetBlockCoordinatesOfFirstIntersectionBlock(iCluster, fromLowerClusterSiteToFirstRayIntersection);

            //The Cluster Traverser keeps track of which block we're at
            //in the cluster
            ClusterTraverser<ClusterType> clusterTraverser(iCluster);

            clusterTraverser.SetCurrentLocation(blockHoldingFirstIntersection);

            // The number of ray units in a block in each direction are cached.
            const util::Vector3D<float> rayUnitsAlongEachBlockSize = ioRay.GetInverseDirection() * blockSizeAsFloat;

            //For every block that is traversed, a vector is needed from
            //where the ray first hits the cluster to the lower site ie site
            //(0,0,0) within the block. This is required to locate how
            //far the ray has travelled and where each part is in relation to
            //voxel sites
            util::Vector3D<float> fromFirstIntersectionToLowerSiteOfCurrentBlock =
                util::Vector3D<float>(blockHoldingFirstIntersection) * blockSizeAsFloat
                    - fromLowerClusterSiteToFirstRayIntersection;

            //We need to know how many ray units can be traversed before
            //a new block is hit. The initial value is calculated based on
            //the location of first intersection
            util::Vector3D<float> totalRayUnitsToNextBlockFromFirstIntersection =
                CalculateMinimalTotalRayUnitsToBlocksBehindCurrentOne(fromFirstIntersectionToLowerSiteOfCurrentBlock,
                                                                      ioRay);

            // We need to track how far the ray has travelled
            float siteUnitsTraversed = 0.0F;

            while (clusterTraverser.CurrentLocationValid())
            {
              // The location of the ray within the block.
              util::Vector3D<float> siteLocationWithinBlock = (ioRay.GetDirection()) * siteUnitsTraversed
                  - fromFirstIntersectionToLowerSiteOfCurrentBlock;

              TraverseRayThroughBlock(fromFirstIntersectionToLowerSiteOfCurrentBlock,
                                      siteLocationWithinBlock,
                                      iCluster,
                                      iCluster.GetMinBlockLocation() + clusterTraverser.GetCurrentLocation(),
                                      clusterTraverser.GetCurrentIndex(),
                                      siteUnitsTraversed,
                                      ioRay);

              // The direction of least travel is the direction of
              // the next block that will be hit by the ray
              // relative to the current block.
              util::Direction::Direction lDirectionOfLeastTravel =
                  DirectionOfLeastTravel(totalRayUnitsToNextBlockFromFirstIntersection);

              //Move to the next block based on the direction
              //of least travel and update variables accordingly
              siteUnitsTraversed =
                  totalRayUnitsToNextBlockFromFirstIntersection.GetByDirection(lDirectionOfLeastTravel);

              switch (lDirectionOfLeastTravel)
              {
                case util::Direction::X:
                  if (ioRay.XIncreasing())
                  {
                    clusterTraverser.IncrementX();
                    fromFirstIntersectionToLowerSiteOfCurrentBlock.x += blockSizeAsFloat;
                    totalRayUnitsToNextBlockFromFirstIntersection.x += rayUnitsAlongEachBlockSize.x;
                  }
                  else
                  {
                    clusterTraverser.DecrementX();
                    fromFirstIntersectionToLowerSiteOfCurrentBlock.x -= blockSizeAsFloat;
                    totalRayUnitsToNextBlockFromFirstIntersection.x -= rayUnitsAlongEachBlockSize.x;
                  }
                  break;

                case util::Direction::Y:
                  if (ioRay.YIncreasing())
                  {
                    clusterTraverser.IncrementY();
                    fromFirstIntersectionToLowerSiteOfCurrentBlock.y += blockSizeAsFloat;
                    totalRayUnitsToNextBlockFromFirstIntersection.y += rayUnitsAlongEachBlockSize.y;
                  }
                  else
                  {
                    clusterTraverser.DecrementY();
                    fromFirstIntersectionToLowerSiteOfCurrentBlock.y -= blockSizeAsFloat;
                    totalRayUnitsToNextBlockFromFirstIntersection.y -= rayUnitsAlongEachBlockSize.y;
                  }
                  break;

                case util::Direction::Z:
                  if (ioRay.ZIncreasing())
                  {
                    clusterTraverser.IncrementZ();
                    fromFirstIntersectionToLowerSiteOfCurrentBlock.z += blockSizeAsFloat;
                    totalRayUnitsToNextBlockFromFirstIntersection.z += rayUnitsAlongEachBlockSize.z;
                  }
                  else
                  {
                    clusterTraverser.DecrementZ();
                    fromFirstIntersectionToLowerSiteOfCurrentBlock.z -= blockSizeAsFloat;
                    totalRayUnitsToNextBlockFromFirstIntersection.z -= rayUnitsAlongEachBlockSize.z;
                  }
                  break;
              }

            }

          }

          util::Vector3D<site_t> GetBlockCoordinatesOfFirstIntersectionBlock(const ClusterType& iCluster,
                                                                             const util::Vector3D<float>& iLowerSiteToFirstRayClusterIntersection)
          {
            util::Vector3D<site_t> lBlockCoordinatesOfFirstIntersectionBlock;

            //Perform the truncated division and ensure that the
            //coordinates are valid to allow for numerical errors
            const util::Vector3D<float> exactBlockCoordsOfFirstIntersectingBlock =
                iLowerSiteToFirstRayClusterIntersection * (1.0F / (float) latticeData.GetBlockSize());

            lBlockCoordinatesOfFirstIntersectionBlock.x =
                (site_t) util::NumericalFunctions::enforceBounds<site_t>((site_t) exactBlockCoordsOfFirstIntersectingBlock.x,
                                                                         0,
                                                                         iCluster.GetBlocksX() - 1);

            lBlockCoordinatesOfFirstIntersectionBlock.y =
                (site_t) util::NumericalFunctions::enforceBounds<site_t>((site_t) exactBlockCoordsOfFirstIntersectingBlock.y,
                                                                         0,
                                                                         iCluster.GetBlocksY() - 1);

            lBlockCoordinatesOfFirstIntersectionBlock.z =
                (site_t) util::NumericalFunctions::enforceBounds<site_t>((site_t) exactBlockCoordsOfFirstIntersectingBlock.z,
                                                                         0,
                                                                         iCluster.GetBlocksZ() - 1);

            return lBlockCoordinatesOfFirstIntersectionBlock;
          }

          util::Vector3D<float> CalculateMinimalTotalRayUnitsToBlocksBehindCurrentOne(const util::Vector3D<float>& lFirstIntersectionToBlockLowerSite,
                                                                                      const Ray<RayDataType>& iRay) const
          {
            util::Vector3D<float> lRayUnits;

            //The ray is currently at the first intersection
            //The number of ray units for a given co-ordinate is the
            //distance in sites to the next block divided by
            //the ray unit distance projected in that direction

            if (iRay.GetDirection().x == 0.0F)
            {
              lRayUnits.x = std::numeric_limits<float>::max();
            }
            else
            {
              lRayUnits.x = lFirstIntersectionToBlockLowerSite.x;
              //If the ray is increasing in this coordinate, we want the
              //distance to the next block
              if (iRay.XIncreasing())
              {
                lRayUnits.x += (float) (latticeData.GetBlockSize());
              }
              //Turn this from site units into ray units
              lRayUnits.x *= iRay.GetInverseDirection().x;
            }

            if (iRay.GetDirection().y == 0.0F)
            {
              lRayUnits.y = std::numeric_limits<float>::max();
            }
            else
            {
              lRayUnits.y = lFirstIntersectionToBlockLowerSite.y;
              if (iRay.YIncreasing())
              {
                lRayUnits.y += (float) (latticeData.GetBlockSize());
              }
              lRayUnits.y *= iRay.GetInverseDirection().y;
            }

            if (iRay.GetDirection().z == 0.0F)
            {
              lRayUnits.z = std::numeric_limits<float>::max();
            }
            else
            {
              lRayUnits.z = lFirstIntersectionToBlockLowerSite.z;
              if (iRay.ZIncreasing())
              {
                lRayUnits.z += (float) (latticeData.GetBlockSize());
              }
              lRayUnits.z *= iRay.GetInverseDirection().z;
            }

            return lRayUnits;
          }

          const Viewpoint& viewpoint;
          const Screen& screen;
          const DomainStats& domainStats;
          const VisSettings& visSettings;
          const hemelb::geometry::LatticeData& latticeData;
          /**
           * The cache of macroscopic fluid properties at each local fluid site.
           */
          const lb::MacroscopicPropertyCache& propertyCache;

          util::Vector3D<float> fromCameraToBottomLeftPixelOfSubImage;

          util::Vector3D<float> mLowerSiteCordinatesOfClusterRelativeToViewpoint;

          XYCoordinates<int> lowerLeftPixelCoordinatesOfSubImage;
          XYCoordinates<int> upperRightPixelCoordinatesOfSubImage;

          //Vectors from the viewpoint centre
          //to the maximum and minimum site span
          //locations respectively
          util::Vector3D<float> mViewpointCentreToMaxSite;
          util::Vector3D<float> mViewpointCentreToMinSite;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERRAYTRACER_H
