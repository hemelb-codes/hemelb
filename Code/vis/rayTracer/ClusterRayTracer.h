#ifndef HEMELB_VIS_CLUSTERRAYTRACER_H
#define HEMELB_VIS_CLUSTERRAYTRACER_H

//#define NDEBUG;
#include <cmath> 
#include <cassert>
#include <iostream>
#include <limits>

#include "geometry/SiteTraverser.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "vis/DomainStats.h"
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
      namespace Direction
      {
        enum Direction
        {
          X,
          Y,
          Z
        };
      }

      template<typename ClusterType, typename RayDataType>
	class ClusterRayTracer
      {
      public:
      ClusterRayTracer(const Viewpoint& iViewpoint,
		       Screen& iScreen,
		       const DomainStats& iDomainStats,
		       const VisSettings& iVisSettings,
		       const hemelb::geometry::LatticeData& iLatticeData) :
	mViewpoint(iViewpoint), mScreen(iScreen), mDomainStats(iDomainStats), mVisSettings(iVisSettings), mLatticeData(iLatticeData)
	{
	}

	void RenderCluster(const ClusterType& iCluster)
	{
	  mLowerSiteCordinatesOfClusterRelativeToViewpoint = iCluster.minBlock
	    - mViewpoint.GetViewpointLocation();

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

	  CastRaysForEachPixel(iCluster);
	}

      protected:
	void GetRayUnitsFromViewpointToCluster(const Ray<RayDataType> & iRay,
					       float & oMaximumRayUnits,
					       float & oMinimumRayUnits)
	{
	  //(Remember that iRay.mDirection is normalised)
	  float lMaxUnitRaysBasedOnX;
	  float lMinUnitRaysBasedOnX;
	  if (iRay.GetDirection().x > 0.0F)
	  {
	    lMaxUnitRaysBasedOnX = mViewpointCentreToMaxSite.x * iRay.GetInverseDirection().x;

	    lMinUnitRaysBasedOnX = mViewpointCentreToMinSite.x * iRay.GetInverseDirection().x;
	  }
	  else if(iRay.GetDirection().x < 0.0F)
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
	  else if(iRay.GetDirection().y < 0.0F)
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
	  else if(iRay.GetDirection().z < 0.0F)
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
	  oMaximumRayUnits = std::min(std::min(lMaxUnitRaysBasedOnX, lMaxUnitRaysBasedOnY),
				      lMaxUnitRaysBasedOnZ);

	  //Maximum ray units to get us into the cluster
	  //We want the maximum number - since only at this point
	  // is the ray completely in
	  oMinimumRayUnits = std::max(std::max(lMinUnitRaysBasedOnX, lMinUnitRaysBasedOnY),
				      lMinUnitRaysBasedOnZ);

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

	  util::Vector3D<float> lLowerSiteToFirstRayClusterIntersection = iMinimumRayUnits
	    * ioRay.GetDirection() - mLowerSiteCordinatesOfClusterRelativeToViewpoint;

	  TraverseBlocks(iCluster, lLowerSiteToFirstRayClusterIntersection, ioRay);
	}

	const Viewpoint& mViewpoint;

	Screen& mScreen;

	const DomainStats& mDomainStats;

	const VisSettings& mVisSettings;

      private:
	void CalculateSubImage(const ClusterType& iCluster)
	{
	  //The extent of the cluster when projected (ie the subimage)
	  //is determined by projecting all eight verticies of the cuboid
	  XYCoordinates<float> lSubImageLowerLeft = XYCoordinates<float>::MaxLimit();
	  XYCoordinates<float> lSubImageUpperRight = XYCoordinates<float>::MinLimit();

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.minSite.x,
							iCluster.minSite.y,
							iCluster.minSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.minSite.x,
							iCluster.minSite.y,
							iCluster.maxSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.minSite.x,
							iCluster.maxSite.y,
							iCluster.minSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.minSite.x,
							iCluster.maxSite.y,
							iCluster.maxSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.maxSite.x,
							iCluster.minSite.y,
							iCluster.minSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.maxSite.x,
							iCluster.minSite.y,
							iCluster.maxSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.maxSite.x,
							iCluster.maxSite.y,
							iCluster.minSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  UpdateSubImageExtentForCorner(util::Vector3D<float>(iCluster.maxSite.x,
							iCluster.maxSite.y,
							iCluster.maxSite.z),
					lSubImageLowerLeft,
					lSubImageUpperRight);

	  mSubImageLowerLeftPixelCoordinates = 
	    mScreen.template TransformScreenToPixelCoordinates<int>(lSubImageLowerLeft);

	  mSubImageUpperRightPixelCoordinates = 
	    mScreen.template TransformScreenToPixelCoordinates<int>(lSubImageUpperRight);
	}

	void UpdateSubImageExtentForCorner(const util::Vector3D<float>& iCorner,
					   XYCoordinates<float>& ioSubImageLowerLeft,
					   XYCoordinates<float>& ioSubImageUpperRight)
	{
	  XYCoordinates<float> lCornerProjection = mViewpoint.FlatProject(iCorner);

	  ioSubImageLowerLeft.UpdatePointwiseMin(lCornerProjection);
	  ioSubImageUpperRight.UpdatePointwiseMax(lCornerProjection);
	}

	bool SubImageOffScreen()
	{
	  if (mSubImageLowerLeftPixelCoordinates.x >= mScreen.GetPixelsX()
	      || mSubImageUpperRightPixelCoordinates.x < 0
	      || mSubImageLowerLeftPixelCoordinates.y >= mScreen.GetPixelsY()
	      || mSubImageUpperRightPixelCoordinates.y < 0)
	  {
	    return true;
	  }
	  else
	  {
	    return false;
	  }
	}

	void CropSubImageToScreen()
	{
	  mSubImageLowerLeftPixelCoordinates.x =
	    util::NumericalFunctions::max(mSubImageLowerLeftPixelCoordinates.x, 0);

	  mSubImageUpperRightPixelCoordinates.x =
	    util::NumericalFunctions::min(mSubImageUpperRightPixelCoordinates.x,
					  mScreen.GetPixelsX() - 1);

	  mSubImageLowerLeftPixelCoordinates.y =
	    util::NumericalFunctions::max(mSubImageLowerLeftPixelCoordinates.y, 0);

	  mSubImageUpperRightPixelCoordinates.y =
	    util::NumericalFunctions::min(mSubImageUpperRightPixelCoordinates.y,
					  mScreen.GetPixelsY() - 1);
	}

	void CalculateVectorsToClusterSpanAndLowerLeftPixel(const ClusterType& iCluster)
	{
	  //This used to be AABB
	  mViewpointCentreToMaxSite = iCluster.maxSite - mViewpoint.GetViewpointLocation();

	  mViewpointCentreToMinSite = iCluster.minSite - mViewpoint.GetViewpointLocation();

	  //Obtaining the vector from the camera to the lower left pixel
	  mCameraToBottomLeftPixel = mScreen.GetCameraToBottomLeftOfScreenVector()
	    + mScreen.GetPixelUnitVectorProjectionX()
	    * (float) mSubImageLowerLeftPixelCoordinates.x
	    + mScreen.GetPixelUnitVectorProjectionY()
	    * (float) mSubImageLowerLeftPixelCoordinates.y;
	}

	void CastRaysForEachPixel(const ClusterType& iCluster)
	{
	  XYCoordinates<int> lPixel;

	  //Loop over all the pixels
	  util::Vector3D<float> lCameraToBottomRow = mCameraToBottomLeftPixel;
	  for (lPixel.x = mSubImageLowerLeftPixelCoordinates.x;
	       lPixel.x <= mSubImageUpperRightPixelCoordinates.x; ++lPixel.x)
	  {
	    util::Vector3D<float> lCameraToPixel = lCameraToBottomRow;
	    for (lPixel.y = mSubImageLowerLeftPixelCoordinates.y;
		 lPixel.y <= mSubImageUpperRightPixelCoordinates.y; ++lPixel.y)
	    {
	      CastRayForPixel(iCluster, lPixel, lCameraToPixel);

	      lCameraToPixel += mScreen.GetPixelUnitVectorProjectionY();
	    }

	    lCameraToBottomRow += mScreen.GetPixelUnitVectorProjectionX();
	  }
	}

	virtual void CastRayForPixel(const ClusterType& iCluster,
				     const XYCoordinates<int>& iPixelCoordinates,
				     const util::Vector3D<float>& iRayDirection)
	{
	  Ray<RayDataType> lRay(iRayDirection);

	  //These tell us how many ray units get us into the cluster
	  //and after how many ray units we are out
	  float lMaximumRayUnits;
	  float lMinimumRayUnits;
	  GetRayUnitsFromViewpointToCluster(lRay, lMaximumRayUnits, lMinimumRayUnits);

	  CastRay(iCluster, lRay, lMaximumRayUnits, lMinimumRayUnits);

	  //Make sure the ray hasn't reached infinity
	  if (!lRay.CollectedNoData())
	  {
	    mScreen.AddRayData(iPixelCoordinates, lRay.GetRayData(), mVisSettings);              
	  }
	}

	void TraverseVoxels(const util::Vector3D<float>& iFirstRayClusterIntersectionToBlockLowerSite,
			    const util::Vector3D<float>& iLocationInBlock,
			    const ClusterType& iCluster,
			    site_t iBlockNumber,
			    float iRayLengthTraversedSoFar,
			    Ray<RayDataType>& ioRay)
	{
	  //Work out which site we're currently in
	  util::Vector3D<site_t> lTruncatedLocationInBlock = RoundToNearestVoxel(iLocationInBlock);

	  geometry::SiteTraverser lSiteTraverser(mLatticeData);
	  lSiteTraverser.SetCurrentLocation(lTruncatedLocationInBlock);

	  //In order to trace the rays through the voxels, we need to
	  //keep track of how far the ray can travel to the next
	  //voxel in each of the three directions in ray units
	  util::Vector3D<float> lRayUnitsBeforeNextVoxel =
	    CalculateRayUnitsBeforeNextVoxel(iFirstRayClusterIntersectionToBlockLowerSite,
					     util::Vector3D<float>(lTruncatedLocationInBlock),
					     ioRay);

	  while (lSiteTraverser.CurrentLocationValid())
	  {
	    //Firstly, work out in which direction we
	    //can travel the least ray units before reaching
	    //a vortex side
	    Direction::Direction lDirectionOfLeastTravel =
	      DirectionOfLeastTravel(lRayUnitsBeforeNextVoxel);

	    float lMinRayUnitsBeforeNextVoxel;
	    //Find out how far the ray can move
	    switch (lDirectionOfLeastTravel)
	    {
	    case Direction::X:
	      lMinRayUnitsBeforeNextVoxel = lRayUnitsBeforeNextVoxel.x;
	      break;

	    case Direction::Y:
	      lMinRayUnitsBeforeNextVoxel = lRayUnitsBeforeNextVoxel.y;
	      break;

	    case Direction::Z:
	      lMinRayUnitsBeforeNextVoxel = lRayUnitsBeforeNextVoxel.z;
	      break;
	    }

	    // Update the ray data
	    // The ray may have been sitting on the fence, in which case
	    // there is little sence in updating ray data
	    if (lMinRayUnitsBeforeNextVoxel - iRayLengthTraversedSoFar > 0.001F)
	    {
	      const SiteData_t* lSiteData = 
		iCluster.GetSiteData(iBlockNumber,
				     lSiteTraverser.GetCurrentIndex());
	      
	      if(lSiteData->Density >= 0.0F) // Ensure fluid site
	      {

		const double* lWallData = iCluster.GetWallData
		  ( iBlockNumber,
		    lSiteTraverser.GetCurrentIndex());
		
		if(lWallData == NULL)
		{
		  ioRay.UpdateDataForNormalFluidSite(*lSiteData,
						     lMinRayUnitsBeforeNextVoxel - iRayLengthTraversedSoFar,
						     iRayLengthTraversedSoFar,
						     mDomainStats,
						     mVisSettings);
		}
		else
		{
		  ioRay.UpdateDataForWallSite(*lSiteData,
					      lMinRayUnitsBeforeNextVoxel - iRayLengthTraversedSoFar,
					      iRayLengthTraversedSoFar,
					      mDomainStats,
					      mVisSettings,
					      lWallData);
		}
	      }
	      
	    }

	    //Update ray length traversed so far
	    iRayLengthTraversedSoFar = lMinRayUnitsBeforeNextVoxel;

	    //Update the block location and RayUnitsBeforeNextVoxel
	    //in each direction
	    switch (lDirectionOfLeastTravel)
	    {
	    case Direction::X:
	      if (ioRay.XIncreasing())
	      {
		lSiteTraverser.IncrementX();
		lRayUnitsBeforeNextVoxel.x += ioRay.GetInverseDirection().x;
	      }
	      else
	      {
		lSiteTraverser.DecrementX();
		lRayUnitsBeforeNextVoxel.x -= ioRay.GetInverseDirection().x;
	      }

	      break;

	    case Direction::Y:
	      if (ioRay.YIncreasing())
	      {
		lSiteTraverser.IncrementY();
		lRayUnitsBeforeNextVoxel.y += ioRay.GetInverseDirection().y;
	      }
	      else
	      {
		lSiteTraverser.DecrementY();
		lRayUnitsBeforeNextVoxel.y -= ioRay.GetInverseDirection().y;
	      }
	      break;

	    case Direction::Z:
	      if (ioRay.ZIncreasing())
	      {
		lSiteTraverser.IncrementZ();
		lRayUnitsBeforeNextVoxel.z += ioRay.GetInverseDirection().z;
	      }
	      else
	      {
		lSiteTraverser.DecrementZ();
		lRayUnitsBeforeNextVoxel.z -= ioRay.GetInverseDirection().z;
	      }
	      break;
	    }

	  }
	}

	util::Vector3D<float> CalculateRayUnitsBeforeNextVoxel
	  (const util::Vector3D<float>& iFirstRayClusterIntersectionToBlockLowerSite,
	   const util::Vector3D<site_t>& iTruncatedLocationInBlock,
	   const Ray<RayDataType>& iRay)
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
	    lRayUnits.x = iFirstRayClusterIntersectionToBlockLowerSite.x
	      + static_cast<float>(iTruncatedLocationInBlock.x);
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
            lRayUnits.y = iFirstRayClusterIntersectionToBlockLowerSite.y
	      + static_cast<float>(iTruncatedLocationInBlock.y);
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
            lRayUnits.z = iFirstRayClusterIntersectionToBlockLowerSite.z
	      + static_cast<float>(iTruncatedLocationInBlock.z);
            if (iRay.ZIncreasing())
            {
              lRayUnits.z += 1.0F;
            }
            lRayUnits.z *= iRay.GetInverseDirection().z;
	  }


	  return lRayUnits;
	}

	util::Vector3D<site_t> RoundToNearestVoxel(const util::Vector3D<float>& iUnboundLocation)
	{
	  util::Vector3D<site_t> lVoxelLocationInBlock;

	  //Due to rounding errors, it's possible for the site location within a block
	  //to be outside the wrong block
	  lVoxelLocationInBlock.x =
	    util::NumericalFunctions::enforceBounds<site_t>((site_t) iUnboundLocation.x,
							    0,
							    mLatticeData.GetBlockSize() - 1);

	  lVoxelLocationInBlock.y =
	    util::NumericalFunctions::enforceBounds<site_t>((site_t) iUnboundLocation.y,
							    0,
							    mLatticeData.GetBlockSize() - 1);

	  lVoxelLocationInBlock.z =
	    util::NumericalFunctions::enforceBounds<site_t>((site_t) iUnboundLocation.z,
							    0,
							    mLatticeData.GetBlockSize() - 1);

	  return lVoxelLocationInBlock;
	}

	Direction::Direction DirectionOfLeastTravel(util::Vector3D<float> iRayUnitsBeforeNextVoxelOrBlock)
	{
	    
	  if (iRayUnitsBeforeNextVoxelOrBlock.x < iRayUnitsBeforeNextVoxelOrBlock.y)
	  {
	    //X is less than Y
	    if (iRayUnitsBeforeNextVoxelOrBlock.x < iRayUnitsBeforeNextVoxelOrBlock.z)
	    {
	      //X is less than Y and X
	      return Direction::X;
	    }
	    else
	    {
	      //X is less than Y
	      //Z is less Than X (And Y)
	      return Direction::Z;
	    }
	  }
	  else
	  {
	    // Y is less than X
	    if (iRayUnitsBeforeNextVoxelOrBlock.y < iRayUnitsBeforeNextVoxelOrBlock.z)
	    {
	      //Y is less than X and Z
	      return Direction::Y;
	    }
	    else
	    {
	      //Y is less than X
	      //Z is less than Y (and so X)
	      return Direction::Z;
	    }

	  }
	}

	void TraverseBlocks(const ClusterType& iCluster,
			    const util::Vector3D<float>& iLowerSiteToFirstRayClusterIntersection,
			    Ray<RayDataType>& ioRay)
	{
	  float lBlockSizeFloat = static_cast<float>(mLatticeData.GetBlockSize());

	  //Calculate the coordinates of the block within the cluster where
	  //the ray first intersects
	  util::Vector3D <site_t> lBlockCoordinatesOfFirstIntersectionBlock =
	    GetBlockCoordinatesOfFirstIntersectionBlock(
	      iCluster, iLowerSiteToFirstRayClusterIntersection);

	  //The Cluster Traverser keeps track of which block we're at
	  //in the cluster
	  ClusterTraverser<ClusterType> lClusterTraverser(iCluster);
	  lClusterTraverser.SetCurrentLocation
	    (lBlockCoordinatesOfFirstIntersectionBlock);

	  //For every block that is traversed, a vector is needed from
	  //where the ray first hits the cluster to the lower site ie site
	  //(0,0,0) within the block. This is required to locate how
	  //far the ray has travelled and where each part is in relation to
	  //voxel sites
	  util::Vector3D<float> lFirstIntersectionToBlockLowerSite =
	    util::Vector3D<float>(lBlockCoordinatesOfFirstIntersectionBlock) * lBlockSizeFloat -
	    iLowerSiteToFirstRayClusterIntersection;

	  //The location where the ray hits each block (relative to site coordintes (0,0,0))
	  //must be calculated to allow correct traversal. Initially this the negative of
	  //the vector between the first intersection in the cluster and site (0,0,0)
	  util::Vector3D <float> lSiteLocationWithinBlock = lFirstIntersectionToBlockLowerSite * -1.0F;

	  //We need to know how many ray units can be traversed before
	  //a new block is hit. The initial value is calculated based on
	  //the location of first intersection
	  util::Vector3D<float> lRayUnitsBeforeNextBlock =
	    CalculateRayUnitsBeforeNextBlock
	    (lFirstIntersectionToBlockLowerSite,
	     ioRay);

	  //The number of ray units in a block in each direction
	  //are cached.
	  util::Vector3D<float> lBlockRayUnitIncrement = ioRay.GetInverseDirection()
	    * lBlockSizeFloat;

	  //We need to track how many site units have been traversed
	  //(ray direction normalises to 1)
	  //from the point of first intersection of the cluster
	  float lSiteUnitsTraversed = 0.0F;

	  while (lClusterTraverser.CurrentLocationValid())
	  {
	    Direction::Direction lDirectionOfLeastTravel =
	      DirectionOfLeastTravel(lRayUnitsBeforeNextBlock);

	    if (iCluster.BlockContainsSites(lClusterTraverser.GetCurrentIndex()))
	    {
	      //Recalculate the site location within the block
	      lSiteLocationWithinBlock =
		lSiteUnitsTraversed * ioRay.GetDirection() -
		lFirstIntersectionToBlockLowerSite;

	      TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			     lSiteLocationWithinBlock,
			     iCluster,
			     lClusterTraverser.GetCurrentIndex(),
			     lSiteUnitsTraversed,
			     ioRay);
	    }

	    //Move to another block
	    switch (lDirectionOfLeastTravel)
	    {
	    case Direction::X:
	      lSiteUnitsTraversed = lRayUnitsBeforeNextBlock.x;
	      if (ioRay.XIncreasing())
	      {
		lClusterTraverser.IncrementX();
		lFirstIntersectionToBlockLowerSite.x += lBlockSizeFloat;
		lRayUnitsBeforeNextBlock.x += lBlockRayUnitIncrement.x;
	      }
	      else
	      {
		lClusterTraverser.DecrementX();
		lFirstIntersectionToBlockLowerSite.x -= lBlockSizeFloat;
		lRayUnitsBeforeNextBlock.x -= lBlockRayUnitIncrement.x;
	      }
	      break;

	    case Direction::Y:
	      lSiteUnitsTraversed = lRayUnitsBeforeNextBlock.y;
	      if (ioRay.YIncreasing())
	      {
		lClusterTraverser.IncrementY();
		lFirstIntersectionToBlockLowerSite.y += lBlockSizeFloat;
		lRayUnitsBeforeNextBlock.y += lBlockRayUnitIncrement.y;
	      }
	      else
	      {
		lClusterTraverser.DecrementY();
		lFirstIntersectionToBlockLowerSite.y -= lBlockSizeFloat;
		lRayUnitsBeforeNextBlock.y -= lBlockRayUnitIncrement.y;
	      }
	      break;

	    case Direction::Z:
	      lSiteUnitsTraversed = lRayUnitsBeforeNextBlock.z;
	      if (ioRay.ZIncreasing())
	      {
		lClusterTraverser.IncrementZ();
		lFirstIntersectionToBlockLowerSite.z += lBlockSizeFloat;
		lRayUnitsBeforeNextBlock.z += lBlockRayUnitIncrement.z;
	      }
	      else
	      {
		lClusterTraverser.DecrementZ();
		lFirstIntersectionToBlockLowerSite.z -= lBlockSizeFloat;
		lRayUnitsBeforeNextBlock.z -= lBlockRayUnitIncrement.z;
	      }
	      break;
	    }

	  }

	}

	util::Vector3D<unsigned int> GetBlockCoordinatesOfFirstIntersectionBlock(
	  const ClusterType& iCluster,
	  util::Vector3D<float> iLowerSiteToFirstRayClusterIntersection)
	{
	  util::Vector3D<unsigned int> lBlockCoordinatesOfFirstIntersectionBlock;

	  //Perform the truncated division and ensure that the
	  //coordinates are valid to allow for numerical errors
	  lBlockCoordinatesOfFirstIntersectionBlock.x = (unsigned int)
	    util::NumericalFunctions::enforceBounds(
	      iCluster.blocksX - 1,
	      0,
	      (int) (1.0F / static_cast<float>(mLatticeData.GetBlockSize())
		     * iLowerSiteToFirstRayClusterIntersection.x));

	  lBlockCoordinatesOfFirstIntersectionBlock.y = (unsigned int)
	    util::NumericalFunctions::enforceBounds(
	      iCluster.blocksY - 1,
	      0,
	      (int) (1.0F / static_cast<float>(mLatticeData.GetBlockSize())
		     * iLowerSiteToFirstRayClusterIntersection.y));

	  lBlockCoordinatesOfFirstIntersectionBlock.z = (unsigned int)
	    util::NumericalFunctions::enforceBounds(
	      iCluster.blocksZ - 1,
	      0,
	      (int) (1.0F / static_cast<float>(mLatticeData.GetBlockSize())
		     * iLowerSiteToFirstRayClusterIntersection.z));

	  return lBlockCoordinatesOfFirstIntersectionBlock;
	}

	util::Vector3D<float> CalculateRayUnitsBeforeNextBlock
	  (const util::Vector3D<float>& lFirstIntersectionToBlockLowerSite,
	   const Ray<RayDataType>& iRay)
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
	    if(iRay.XIncreasing())
	    {
	      lRayUnits.x += static_cast<float>(mLatticeData.GetBlockSize());
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
	    if(iRay.YIncreasing())
	    {
	      lRayUnits.y += static_cast<float>(mLatticeData.GetBlockSize());
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
	    if(iRay.ZIncreasing())
	    {
	      lRayUnits.z += static_cast<float>(mLatticeData.GetBlockSize());
	    }
	    lRayUnits.z *= iRay.GetInverseDirection().z;
	  }
		
	  return lRayUnits;
	}

        
	const hemelb::geometry::LatticeData& mLatticeData;

	util::Vector3D<float> mCameraToBottomLeftPixel;

	util::Vector3D<float> mLowerSiteCordinatesOfClusterRelativeToViewpoint;

	XYCoordinates<int> mSubImageLowerLeftPixelCoordinates;

	XYCoordinates<int> mSubImageUpperRightPixelCoordinates;

	//Vectors from the viewpoint centre
	//to the maximum and minimum site span
	//locations respectively
	//(Formerly AABB)
	util::Vector3D<float> mViewpointCentreToMaxSite;
	util::Vector3D<float> mViewpointCentreToMinSite;
      };}
  }
}

#endif // HEMELB_VIS_CLUSTERRENDERER_H
