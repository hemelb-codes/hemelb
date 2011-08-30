//#define NDEBUG;
#include <cmath> 
#include <cassert>
#include <iostream>

#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/ClusterTraverser.h"
#include "vis/rayTracer/ClusterRayTracer.h"
#include "vis/rayTracer/Ray.h"
#include "vis/rayTracer/SiteTraverser.h"

#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      ClusterRayTracer::ClusterRayTracer
      (	const Viewpoint& iViewpoint,
	Screen& iScreen,
	const DomainStats& iDomainStats,
	const VisSettings& iVisSettings,
	const hemelb::geometry::LatticeData& iLatticeData) :
	mViewpoint(iViewpoint),
	mScreen(iScreen),
	mDomainStats(iDomainStats),
	mVisSettings(iVisSettings),
	mLatticeData(iLatticeData)
      {
      }

      void ClusterRayTracer::RenderCluster(const Cluster& iCluster)
      {
	mLowerSiteCordinatesOfClusterRelativeToViewpoint = 
	  iCluster.minBlock - mViewpoint.GetViewpointLocation();
	
	//Calculate the projection of the cluster on the screen
	//refered to as the subimage
	CalculateSubImage(iCluster);
	
	//If the entire sub-image is off the screen,
        //no rendering is needed
	if(SubImageOffScreen())
	{
	  return;
	}

	CropSubImageToScreen();

	CalculateVectorsToClusterSpanAndLowerLeftPixel(iCluster);

	CastRaysForEachPixel(iCluster);
      }

      void ClusterRayTracer::CalculateSubImage(const Cluster& iCluster)
      {
	//The extent of the cluster when projected (ie the subimage) 
	//is determined by projecting all eight verticies of the cuboid
	XYCoordinates<float> lSubImageLowerLeft = XYCoordinates<float>::MaxLimit();
	XYCoordinates<float> lSubImageUpperRight = XYCoordinates<float>::MinLimit();

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.minSite.x,
					iCluster.minSite.y,
					iCluster.minSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.minSite.x,
					iCluster.minSite.y,
					iCluster.maxSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.minSite.x,
					iCluster.maxSite.y,
					iCluster.minSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.minSite.x,
					iCluster.maxSite.y,
					iCluster.maxSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.maxSite.x,
					iCluster.minSite.y,
					iCluster.minSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.maxSite.x,
					iCluster.minSite.y,
					iCluster.maxSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.maxSite.x,
					iCluster.maxSite.y,
					iCluster.minSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	UpdateSubImageExtentForCorner(Vector3D<float>(
					iCluster.maxSite.x,
					iCluster.maxSite.y,
					iCluster.maxSite.z),
				      lSubImageLowerLeft,
				      lSubImageUpperRight);

	mSubImageLowerLeftPixelCoordinates = 
	  mScreen.TransformScreenToPixelCoordinates<int> (lSubImageLowerLeft);

	mSubImageUpperRightPixelCoordinates = 
	  mScreen.TransformScreenToPixelCoordinates<int> (lSubImageUpperRight);
	
      }

      void ClusterRayTracer::UpdateSubImageExtentForCorner
      (const Vector3D<float>& iCorner,
       XYCoordinates<float>& ioSubImageLowerLeft,
       XYCoordinates<float>& ioSubImageUpperRight)
      {
	XYCoordinates<float> lCornerProjection = mViewpoint.FlatProject(iCorner);
	
	XYCoordinates<float>::UpdateMinXYCoordinates(ioSubImageLowerLeft, lCornerProjection);
	XYCoordinates<float>::UpdateMaxXYCoordinates(ioSubImageUpperRight, lCornerProjection);
      }

      bool ClusterRayTracer::SubImageOffScreen()
      {
	if (mSubImageLowerLeftPixelCoordinates.x >= mScreen.GetPixelsX() ||
	    mSubImageUpperRightPixelCoordinates.x < 0 ||
	    mSubImageLowerLeftPixelCoordinates.y >= mScreen.GetPixelsY() ||
	    mSubImageUpperRightPixelCoordinates.y < 0)
	{
	   std::cout << "Off screen" << std::endl;
	   return true;
	}
	else
	{
	  return false;
	}
      }

      void ClusterRayTracer::CropSubImageToScreen()
      {
	mSubImageLowerLeftPixelCoordinates.x = util::NumericalFunctions::
	  max(mSubImageLowerLeftPixelCoordinates.x, 0);
	
	mSubImageUpperRightPixelCoordinates.x = util::NumericalFunctions::
	  min(mSubImageUpperRightPixelCoordinates.x, mScreen.GetPixelsX() - 1);
	
	mSubImageLowerLeftPixelCoordinates.y = util::NumericalFunctions::
	  max(mSubImageLowerLeftPixelCoordinates.y, 0);

	mSubImageUpperRightPixelCoordinates.y = util::NumericalFunctions::
	  min(mSubImageUpperRightPixelCoordinates.y, mScreen.GetPixelsY() - 1);
      }

      void ClusterRayTracer::
      CalculateVectorsToClusterSpanAndLowerLeftPixel(const Cluster& iCluster)
      {
	//This used to be AABB
	mViewpointCentreToMaxSite = iCluster.maxSite - 
	  mViewpoint.GetViewpointLocation();

	mViewpointCentreToMinSite = iCluster.minSite - 
	  mViewpoint.GetViewpointLocation();

	//Obtaining the vector from the camera to the lower left pixel
	mCameraToBottomLeftPixel = mScreen.GetCameraToBottomLeftOfScreenVector() + 
	  mScreen.GetPixelUnitVectorProjectionX() * (float) mSubImageLowerLeftPixelCoordinates.x +
	  mScreen.GetPixelUnitVectorProjectionY() * (float) mSubImageLowerLeftPixelCoordinates.y;
      }

      void ClusterRayTracer::CastRaysForEachPixel(const Cluster& iCluster)
      {
	XYCoordinates<int> lPixel;

	//Loop over all the pixels
	Vector3D<float> lCameraToBottomRow = mCameraToBottomLeftPixel;
	for (lPixel.x = mSubImageLowerLeftPixelCoordinates.x;
	     lPixel.x <= mSubImageUpperRightPixelCoordinates.x;
	     ++lPixel.x)
	{
	  Vector3D<float> lCameraToPixel = lCameraToBottomRow; 
	  for (lPixel.y = mSubImageLowerLeftPixelCoordinates.y;
	       lPixel.y <= mSubImageUpperRightPixelCoordinates.y;
	       ++lPixel.y)
	  {
	    CastRayForPixel(iCluster, lPixel, lCameraToPixel);

	    lCameraToPixel += mScreen.GetPixelUnitVectorProjectionY();
	  }
	  
	  lCameraToBottomRow+= mScreen.GetPixelUnitVectorProjectionX();
	}
      }


      void ClusterRayTracer::CastRayForPixel(const Cluster& iCluster, 
					     const XYCoordinates<int>& iPixel,
					     Vector3D<float> iRayDirection)
      {
	//This changes the vector hence the vector must be a copy
	iRayDirection.Normalise();

	Ray lRay(iRayDirection);

	//These tell us how many ray units get us into the cluster
	//and after how many ray units we are out
	float lMaximumRayUnits;
	float lMinimumRayUnits;
	GetRayUnitsFromViewpointToCluster
	  (lRay, lMaximumRayUnits, lMinimumRayUnits);

	//It's possible for the ray to totally miss the cluster
	//This is because the sub-image is square while the cluster
	// projection won't be in most circumstances
	if(lMaximumRayUnits < lMinimumRayUnits)
	{
	  return;
	}

	Vector3D <float> lLowerSiteToFirstRayClusterIntersection = 
	  lMinimumRayUnits * lRay.GetDirection() -
	  mLowerSiteCordinatesOfClusterRelativeToViewpoint;

	TraverseBlocks(iCluster, lLowerSiteToFirstRayClusterIntersection, lRay);

	if (lRay.MinT == std::numeric_limits<float>::max())
	{
	  return;
	}

	ColPixel col_pixel(iPixel.x, iPixel.y, lRay.MinT + lMinimumRayUnits, lRay.Length, 
			   (lRay.Density - (float) mDomainStats.density_threshold_min)
			   * (float) mDomainStats.density_threshold_minmax_inv, lRay.Stress
			   != std::numeric_limits<float>::max()
			   ? lRay.Stress * (float) mDomainStats.stress_threshold_max_inv
			   : std::numeric_limits<float>::max(), lRay.VelocityColour, lRay.StressColour);

	mScreen.AddPixel(&col_pixel, &mVisSettings);
      }


      void ClusterRayTracer::GetRayUnitsFromViewpointToCluster
      (const Ray & iRay, 
       float & oMaximumRayUnits, float & oMinimumRayUnits)
      {
	//(Remember that iRay.mDirection is normalised)
	float lMaxUnitRaysBasedOnX;
	float lMinUnitRaysBasedOnX;
	if (iRay.GetDirection().x > 0.0F)
	{
	  lMaxUnitRaysBasedOnX = 
	    mViewpointCentreToMaxSite.x * iRay.GetInverseDirection().x;

	  lMinUnitRaysBasedOnX = 
	    mViewpointCentreToMinSite.x * iRay.GetInverseDirection().x;
	}
	else
	{
	  lMaxUnitRaysBasedOnX = 
	    mViewpointCentreToMinSite.x * iRay.GetInverseDirection().x;
	  lMinUnitRaysBasedOnX = 
	    mViewpointCentreToMaxSite.x * iRay.GetInverseDirection().x;
	}

	float lMaxUnitRaysBasedOnY;
	float lMinUnitRaysBasedOnY;
	if (iRay.GetDirection().y > 0.0F)
	{
	  lMaxUnitRaysBasedOnY = 
	    mViewpointCentreToMaxSite.y * iRay.GetInverseDirection().y;

	  lMinUnitRaysBasedOnY = 
	    mViewpointCentreToMinSite.y * iRay.GetInverseDirection().y;
	}
	else
	{
	  lMaxUnitRaysBasedOnY = 
	    mViewpointCentreToMinSite.y * iRay.GetInverseDirection().y;
	  lMinUnitRaysBasedOnY = 
	    mViewpointCentreToMaxSite.y * iRay.GetInverseDirection().y;
	}
	
	float lMaxUnitRaysBasedOnZ;
	float lMinUnitRaysBasedOnZ;
	if (iRay.GetDirection().z > 0.0F)
	{
	  lMaxUnitRaysBasedOnZ = 
	    mViewpointCentreToMaxSite.z * iRay.GetInverseDirection().z;

	  lMinUnitRaysBasedOnZ = 
	    mViewpointCentreToMinSite.z * iRay.GetInverseDirection().z;
	}
	else
	{
	  lMaxUnitRaysBasedOnZ = 
	    mViewpointCentreToMinSite.z * iRay.GetInverseDirection().z;
	  lMinUnitRaysBasedOnZ = 
	    mViewpointCentreToMaxSite.z * iRay.GetInverseDirection().z;
	}

	//Maximum ray units from viewpoint to cluster
	//We want the minimum number - since at this point the ray is
	//completely out
	oMaximumRayUnits = std::min(std::min(lMaxUnitRaysBasedOnX,lMaxUnitRaysBasedOnY),
			       lMaxUnitRaysBasedOnZ);
	
	//Maximum ray units to get us into the cluster
	//We want the maximum number - since only at this point 
	// is the ray completely in
	oMinimumRayUnits = std::max(std::max(lMinUnitRaysBasedOnX,lMinUnitRaysBasedOnY),
			       lMinUnitRaysBasedOnZ);
      }

      void ClusterRayTracer::TraverseVoxels(
	const Vector3D<float>& iFirstRayClusterIntersectionToBlockLowerSite,
	const Vector3D<float>& iLocationInBlock,
	const SiteData_t* iSiteData,
	float iRayLengthTraversedSoFar,
	Ray& iRay)
      {
	//Work out which site we're currently in
	Vector3D<site_t> lTruncatedLocationInBlock = RoundToNearestVoxel(iLocationInBlock);
	
	SiteTraverser lSiteTraverser(&mLatticeData);
	lSiteTraverser.SetCurrentLocation(lTruncatedLocationInBlock);

	//In order to trace the rays through the voxels, we need to 
	//keep track of how far the ray can travel to the next
	//voxel in each of the three directions in ray units
	Vector3D<float> lRayUnitsBeforeNextVoxel =
	  CalculateRayUnitsBeforeNextVoxel(iFirstRayClusterIntersectionToBlockLowerSite, 
					   Vector3D<float>(lTruncatedLocationInBlock),
					   iRay);

	while(lSiteTraverser.CurrentLocationValid())
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
	  if(lMinRayUnitsBeforeNextVoxel - iRayLengthTraversedSoFar > 0.001F)
	  {
	    UpdateRayData(&iSiteData[lSiteTraverser.GetCurrentIndex()],
			  iRayLengthTraversedSoFar,
			  lMinRayUnitsBeforeNextVoxel - iRayLengthTraversedSoFar,
			  &iRay);
	    }
	  
	  //Update ray length traversed so far
	  iRayLengthTraversedSoFar = lMinRayUnitsBeforeNextVoxel;

	  //Update the block location and RayUnitsBeforeNextVoxel
	  //in each direction
	  switch (lDirectionOfLeastTravel)
	  {
	  case Direction::X:
	    if (iRay.XIncreasing())
	    {
	      lSiteTraverser.IncrementX();
	      lRayUnitsBeforeNextVoxel.x += iRay.GetInverseDirection().x;
	    }
	    else
	    {
	      lSiteTraverser.DecrementX();
	      lRayUnitsBeforeNextVoxel.x -= iRay.GetInverseDirection().x;
	    }

	    break;

	  case Direction::Y:
	    if (iRay.YIncreasing())
	    {
	      lSiteTraverser.IncrementY();
	      lRayUnitsBeforeNextVoxel.y += iRay.GetInverseDirection().y;
	    }
	    else
	    {
	      lSiteTraverser.DecrementY();
	      lRayUnitsBeforeNextVoxel.y -= iRay.GetInverseDirection().y;
	    }
	    break;
	  
	  case Direction::Z:
	    if (iRay.ZIncreasing())
	    {
	      lSiteTraverser.IncrementZ();
	      lRayUnitsBeforeNextVoxel.z += iRay.GetInverseDirection().z;
	    }
	    else
	    {
	      lSiteTraverser.DecrementZ();
	      lRayUnitsBeforeNextVoxel.z -= iRay.GetInverseDirection().z;
	    }   
	    break;
	  }

	}
      }

      Vector3D<site_t> ClusterRayTracer::
      RoundToNearestVoxel(const Vector3D<float>& iUnboundLocation)
      {
	Vector3D<site_t> lVoxelLocationInBlock;
	
	//Due to rounding errors, it's possible for the site location within a block
	//to be outside the wrong block
	lVoxelLocationInBlock.x = util::NumericalFunctions::
	enforceBounds<site_t>((site_t) iUnboundLocation.x,
				0,
				mLatticeData.GetBlockSize() - 1);
	  
	lVoxelLocationInBlock.y  = util::NumericalFunctions::
	  enforceBounds<site_t>((site_t) iUnboundLocation.y,
				0,
				mLatticeData.GetBlockSize() - 1);

	lVoxelLocationInBlock.z = util::NumericalFunctions::
	  enforceBounds<site_t>((site_t) iUnboundLocation.z,
				0,
				mLatticeData.GetBlockSize() - 1);

	return lVoxelLocationInBlock;
      }
 
      Vector3D<float> ClusterRayTracer::CalculateRayUnitsBeforeNextVoxel
      (const Vector3D<float>& iFirstRayClusterIntersectionToBlockLowerSite,
       const Vector3D<site_t>& iTruncatedLocationInBlock, const Ray& iRay)
      {
	Vector3D<float> lRayUnits;
	
	//The ray has already travelled iFirstRayClusterIntersectionToBlockLowerSite
        //If the ray is increasing in the co-ordinte it can travel as far
	//as the truncated location + 1, otherwise just the truncated location
        //for each co-ordinate

	lRayUnits.x = iFirstRayClusterIntersectionToBlockLowerSite.x
	  + static_cast<float>(iTruncatedLocationInBlock.x);
	if(iRay.XIncreasing())
	{
	  lRayUnits.x += 1.0F;
	}
	//Convert from site units into ray units
	lRayUnits.x *= iRay.GetInverseDirection().x;

	lRayUnits.y = iFirstRayClusterIntersectionToBlockLowerSite.y
	  + static_cast<float>(iTruncatedLocationInBlock.y);
	if(iRay.YIncreasing())
	{
	  lRayUnits.y += 1.0F;
	}
	lRayUnits.y *= iRay.GetInverseDirection().y;

	lRayUnits.z = iFirstRayClusterIntersectionToBlockLowerSite.z
	  + static_cast<float>(iTruncatedLocationInBlock.z);
	if(iRay.ZIncreasing())
	{
	  lRayUnits.z += 1.0F;
	}
	lRayUnits.z *= iRay.GetInverseDirection().z;

	return lRayUnits;
      }

    Direction::Direction ClusterRayTracer::
    DirectionOfLeastTravel
    (Vector3D<float> iRayUnitsBeforeNextVoxelOrBlock)
    {
      if (iRayUnitsBeforeNextVoxelOrBlock.x <
	iRayUnitsBeforeNextVoxelOrBlock.y)
      {
	//X is less than Y
	if (iRayUnitsBeforeNextVoxelOrBlock.x <
	iRayUnitsBeforeNextVoxelOrBlock.z)
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
	if (iRayUnitsBeforeNextVoxelOrBlock.y <
	iRayUnitsBeforeNextVoxelOrBlock.z)
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
    

      void ClusterRayTracer::TraverseBlocks(const Cluster& iCluster, 
					    const Vector3D<float>& iLowerSiteToFirstRayClusterIntersection,
					    Ray& iRay)
      {
	float lBlockSizeFloat = static_cast<float>(mLatticeData.GetBlockSize());
	
	//Calculate the coordinates of the block within the cluster where
	//the ray first intersects
 	Vector3D<site_t> lBlockCoordinatesOfFirstIntersectionBlock =  
	  GetBlockCoordinatesOfFirstIntersectionBlock(
	    iCluster, iLowerSiteToFirstRayClusterIntersection);
	
	//The Cluster Traverser keeps track of which block we're at
	//in the cluster
	ClusterTraverser lClusterTraverser(iCluster);
	lClusterTraverser.SetCurrentLocation
	  (lBlockCoordinatesOfFirstIntersectionBlock);

	//For every block that is traversed, a vector is needed from
	//where the ray first hits the cluster to the lower site ie site
	//(0,0,0) within the block. This is required to locate how
	//far the ray has travelled and where each part is in relation to
	//voxel sites
	Vector3D<float> lFirstIntersectionToBlockLowerSite = 
	  Vector3D<float>(lBlockCoordinatesOfFirstIntersectionBlock) * lBlockSizeFloat - 
	  iLowerSiteToFirstRayClusterIntersection;

	//The location where the ray hits each block (relative to site coordintes (0,0,0))
	//must be calculated to allow correct traversal. Initially this the negative of 
	//the vector between the first intersection in the cluster and site (0,0,0)
       	Vector3D <float> lSiteLocationWithinBlock = lFirstIntersectionToBlockLowerSite * -1.0F;
	
	//We need to know how many ray units can be traversed before
	//a new blockis hit. The initial value is calculated based on
	//the location of first intersection
	Vector3D<float> lRayUnitsBeforeNextBlock = 
	  CalculateRayUnitsBeforeNextBlock
	  (lFirstIntersectionToBlockLowerSite,
	   iRay);

	//The number of ray units in a block in each direction 
	//are cached.
	Vector3D<float> lBlockRayUnitIncrement = iRay.GetInverseDirection() 
	  * lBlockSizeFloat;
	
	//We need to track how many ray units have been traversed 
	//from the point of first intersection
	float lRayUnitsTraversed = 0.0F;

	while (lClusterTraverser.CurrentLocationValid())
	{
	  Direction::Direction lDirectionOfLeastTravel =
	    DirectionOfLeastTravel(lRayUnitsBeforeNextBlock);

	  if (iCluster.BlockContainsSites(lClusterTraverser.GetCurrentIndex()))
	  {
	    //Recalculate the site location within the block
	    lSiteLocationWithinBlock = 
	    lRayUnitsTraversed * iRay.GetDirection() -
	    lFirstIntersectionToBlockLowerSite;

	    TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			   lSiteLocationWithinBlock,
			   iCluster.GetSiteData
			   (lClusterTraverser.GetCurrentIndex()),
			   lRayUnitsTraversed,
			   iRay);
	  }

	  //Move to another block
	  switch (lDirectionOfLeastTravel)
	  {
	  case Direction::X:
	    lRayUnitsTraversed = lRayUnitsBeforeNextBlock.x;
	    if (iRay.XIncreasing())
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
	    lRayUnitsTraversed = lRayUnitsBeforeNextBlock.y;
	    if (iRay.YIncreasing())
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
	    lRayUnitsTraversed = lRayUnitsBeforeNextBlock.z;
	    if (iRay.ZIncreasing())
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

      
      Vector3D<unsigned int>  ClusterRayTracer::GetBlockCoordinatesOfFirstIntersectionBlock(
	const Cluster& iCluster,
	Vector3D<float> iLowerSiteToFirstRayClusterIntersection)
      {
	Vector3D<unsigned int> lBlockCoordinatesOfFirstIntersectionBlock;

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
      
      Vector3D<float> ClusterRayTracer::CalculateRayUnitsBeforeNextBlock
      (const Vector3D<float>& lFirstIntersectionToBlockLowerSite,
       const Ray& iRay)
      {
	Vector3D<float> lRayUnits;
	
	//The ray is currently at the first intersection
	//The number of ray units for a given co-ordinate is the
	//distance in sites to the next block divided by
	//the ray unit distance projected in that direction
	
	lRayUnits.x = lFirstIntersectionToBlockLowerSite.x;
	//If the ray is increasing in this coordinate, we want the
	//distance to the next block
	if(iRay.XIncreasing())
	{
	  lRayUnits.x += static_cast<float>(mLatticeData.GetBlockSize());
	}
	//Turn this from site units into ray units
	lRayUnits.x *= iRay.GetInverseDirection().x;
	
	lRayUnits.y = lFirstIntersectionToBlockLowerSite.y;
	if(iRay.YIncreasing())
	{
	  lRayUnits.y += static_cast<float>(mLatticeData.GetBlockSize());
	}
	lRayUnits.y *= iRay.GetInverseDirection().y;


	lRayUnits.z = lFirstIntersectionToBlockLowerSite.z;
	if(iRay.ZIncreasing())
	{
	  lRayUnits.z += static_cast<float>(mLatticeData.GetBlockSize());
	}
	lRayUnits.z *= iRay.GetInverseDirection().z;

	return lRayUnits;
      }


      void ClusterRayTracer::UpdateRayData(const SiteData_t* iSiteData,
				    float ray_t,
				    float ray_segment,
				    Ray* iRay)
      {
	if (iSiteData->Density < 0.0F)
	{
	  return; // solid voxel
	}

	float palette[3];

	// update the volume rendering of the velocity flow field
	ColPixel::PickColour(iSiteData->Velocity * (float) mDomainStats.velocity_threshold_max_inv,
			     palette);

	UpdateColour(ray_segment, palette, iRay->VelocityColour);

	if (mVisSettings.mStressType != lb::ShearStress)
	{
	  // update the volume rendering of the von Mises stress flow field
	  float scaled_stress = iSiteData->Stress * (float) mDomainStats.stress_threshold_max_inv;

	  ColPixel::PickColour(scaled_stress, palette);

	  UpdateColour(ray_segment, palette, iRay->StressColour);
	}

	iRay->Length += ray_segment;

	if (iRay->Density >= 0.0F)
	{
	  return;
	}

	iRay->MinT = ray_t;

	// keep track of the density nearest to the view point
	iRay->Density = iSiteData->Density;

	// keep track of the stress nearest to the view point
	iRay->Stress = iSiteData->Stress;		
      }

      /**
       * Update a colour vector to include a section of known length through a
       * solid of known colour.
       *
       * @param dt
       * @param palette
       * @param col
       */
      void ClusterRayTracer::UpdateColour(float dt, const float palette[3], float col[3])
      {
	for (int ii = 0; ii < 3; ++ii)
	{
	  col[ii] += dt * palette[ii];
	}
      }

    }
  }
}
