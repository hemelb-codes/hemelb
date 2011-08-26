//#define NDEBUG;
#include <cmath> 
#include <cassert>
#include <iostream>

#include "vis/rayTracer/Cluster.h"
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


	Vector3D<bool> lRayInPositiveDirection;
	lRayInPositiveDirection.x = lRay.GetDirection().x > 0.0F;
	lRayInPositiveDirection.y = lRay.GetDirection().y > 0.0F;
	lRayInPositiveDirection.z = lRay.GetDirection().z > 0.0F; 

	TraverseBlocks(iCluster, lRayInPositiveDirection, lLowerSiteToFirstRayClusterIntersection, &lRay);

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
	Ray* bCurrentRay)
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
					   *bCurrentRay);

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
			  bCurrentRay);
	    }
	  
	  //Update ray length traversed so far
	  iRayLengthTraversedSoFar = lMinRayUnitsBeforeNextVoxel;

	  //Update the block location and RayUnitsBeforeNextVoxel
	  //in each direction
	  switch (lDirectionOfLeastTravel)
	  {
	  case Direction::X:
	    if (bCurrentRay->XIncreasing())
	    {
	      lSiteTraverser.IncrementX();
	      lRayUnitsBeforeNextVoxel.x += bCurrentRay->GetInverseDirection().x;
	    }
	    else
	    {
	      lSiteTraverser.DecrementX();
	      lRayUnitsBeforeNextVoxel.x -= bCurrentRay->GetInverseDirection().x;
	    }

	    break;

	  case Direction::Y:
	    if (bCurrentRay->YIncreasing())
	    {
	      lSiteTraverser.IncrementY();
	      lRayUnitsBeforeNextVoxel.y += bCurrentRay->GetInverseDirection().y;
	    }
	    else
	    {
	      lSiteTraverser.DecrementY();
	      lRayUnitsBeforeNextVoxel.y -= bCurrentRay->GetInverseDirection().y;
	    }
	    break;
	  
	  case Direction::Z:
	    if (bCurrentRay->ZIncreasing())
	    {
	      lSiteTraverser.IncrementZ();
	      lRayUnitsBeforeNextVoxel.z += bCurrentRay->GetInverseDirection().z;
	    }
	    else
	    {
	      lSiteTraverser.DecrementZ();
	      lRayUnitsBeforeNextVoxel.z -= bCurrentRay->GetInverseDirection().z;
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
    (Vector3D<float> iRayUnitsBeforeNextVoxel)
    {
      if (iRayUnitsBeforeNextVoxel.x <
	iRayUnitsBeforeNextVoxel.y)
      {
	//X is less than Y
	if (iRayUnitsBeforeNextVoxel.x <
	iRayUnitsBeforeNextVoxel.z)
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
	if (iRayUnitsBeforeNextVoxel.y <
	iRayUnitsBeforeNextVoxel.z)
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
					    const Vector3D<bool>& xyz_Is_1,
					    const Vector3D<float>& iLowerSiteToFirstRayClusterIntersection,
				     Ray *bCurrentRay)
      {
	int cluster_blocksZ = iCluster.blocksZ;
	int cluster_blocksYz = (int) iCluster.blocksY * (int) iCluster.blocksZ;
	int cluster_blocks = (int) iCluster.blocksX * cluster_blocksYz;

	Vector3D<float> iFirstRayClusterIntersectionToBlockLowerSite;

	float lBlockSizeFloat = static_cast<float>(mLatticeData.GetBlockSize());

 	Vector3D<unsigned int> lBlockCoordinatesOfFirstIntersectionBlock = 
	  GetBlockCoordinatesOfFirstIntersectionBlock(iCluster, iLowerSiteToFirstRayClusterIntersection);
	
	Vector3D<float> lFirstIntersectionToBlockLowerSite = 
	  Vector3D<float>(lBlockCoordinatesOfFirstIntersectionBlock) * lBlockSizeFloat - 
	  iLowerSiteToFirstRayClusterIntersection;
	
	int i = lBlockCoordinatesOfFirstIntersectionBlock.x * cluster_blocksYz;
	int j = lBlockCoordinatesOfFirstIntersectionBlock.y * cluster_blocksZ;
	int k = lBlockCoordinatesOfFirstIntersectionBlock.z; 

	//unsigned int lBlockId = iCluster.GetBlockIdFrom3DBlockLocation(lBlockCoordinatesOfFirstIntersectionBlock);

	Vector3D<float> lSiteLocationWithinBlock;
	if (!iCluster.SiteData[i + j + k].empty())
	{
	  lSiteLocationWithinBlock = lFirstIntersectionToBlockLowerSite * -1.0F;

	  TraverseVoxels(lFirstIntersectionToBlockLowerSite, lSiteLocationWithinBlock, &iCluster.SiteData[i + j + k][0], 0.0F, bCurrentRay);
	}

	Vector3D <float> t_max;

	t_max.x = (xyz_Is_1.x
		   ? lFirstIntersectionToBlockLowerSite.x + lBlockSizeFloat
		   : lFirstIntersectionToBlockLowerSite.x) * 1.0F * bCurrentRay->GetInverseDirection().x;

	t_max.y = (xyz_Is_1.y
		   ? lFirstIntersectionToBlockLowerSite.y + lBlockSizeFloat
		   : lFirstIntersectionToBlockLowerSite.y) * 1.0F * bCurrentRay->GetInverseDirection().y;

	t_max.z = (xyz_Is_1.z
		   ? lFirstIntersectionToBlockLowerSite.z + lBlockSizeFloat
		   : lFirstIntersectionToBlockLowerSite.z) * 1.0F * bCurrentRay->GetInverseDirection().z;


	Vector3D<float> t_delta = 1.0F * bCurrentRay->GetInverseDirection() * lBlockSizeFloat;
	
	while (true)
	{
	  if (t_max.x < t_max.y)
	  {
	    if (t_max.x < t_max.z)
	    {
	      if (xyz_Is_1.x)
	      {
		if ( (i += cluster_blocksYz) >= cluster_blocks)
		  return;
		lFirstIntersectionToBlockLowerSite.x += lBlockSizeFloat;
	      }
	      else
	      {
		if ( (i -= cluster_blocksYz) < 0)
		  return;
		lFirstIntersectionToBlockLowerSite.x -= lBlockSizeFloat;
	      }

	      if (!iCluster.SiteData[i + j + k].empty())
	      {
		lSiteLocationWithinBlock.x = t_max.x * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		lSiteLocationWithinBlock.y = t_max.x * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		lSiteLocationWithinBlock.z = t_max.x * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       lSiteLocationWithinBlock,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.x,
			       bCurrentRay);
	      }

	      t_max.x = xyz_Is_1.x
		? t_max.x + t_delta.x
		: t_max.x - t_delta.x;
	    }
	    else
	    {
	      if (xyz_Is_1.z)
	      {
		if (++k >= cluster_blocksZ)
		  return;
		lFirstIntersectionToBlockLowerSite.z += lBlockSizeFloat;
	      }
	      else
	      {
		if (--k < 0)
		  return;
		lFirstIntersectionToBlockLowerSite.z -= lBlockSizeFloat;
	      }

	      if (!iCluster.SiteData[i + j + k].empty())
	      {
		lSiteLocationWithinBlock.x = t_max.z * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		lSiteLocationWithinBlock.y = t_max.z * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		lSiteLocationWithinBlock.z = t_max.z * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       lSiteLocationWithinBlock,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.z,
			       bCurrentRay);
	      }

	      t_max.z = xyz_Is_1.z
		? t_max.z + t_delta.z
		: t_max.z - t_delta.z;
	    }
	  }
	  else
	  {
	    if (t_max.y < t_max.z)
	    {
	      if (xyz_Is_1.y)
	      {
		if ( (j += cluster_blocksZ) >= cluster_blocksYz)
		  return;
		lFirstIntersectionToBlockLowerSite.y += lBlockSizeFloat;
	      }
	      else
	      {
		if ( (j -= cluster_blocksZ) < 0)
		  return;
		lFirstIntersectionToBlockLowerSite.y -= lBlockSizeFloat;
	      }

	      if (!iCluster.SiteData[i + j + k].empty())
	      {
		lSiteLocationWithinBlock.x = t_max.y * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		lSiteLocationWithinBlock.y = t_max.y * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		lSiteLocationWithinBlock.z = t_max.y * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       lSiteLocationWithinBlock,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.y,
			       bCurrentRay);
	      }

	      t_max.y = xyz_Is_1.y
		? t_max.y + t_delta.y
		: t_max.y - t_delta.y;
	    }
	    else
	    {
	      if (xyz_Is_1.z)
	      {
		if (++k >= cluster_blocksZ)
		  return;
		lFirstIntersectionToBlockLowerSite.z += lBlockSizeFloat;
	      }
	      else
	      {
		if (--k < 0)
		  return;
		lFirstIntersectionToBlockLowerSite.z -= lBlockSizeFloat;
	      }

	      if (!iCluster.SiteData[i + j + k].empty())
	      {
		lSiteLocationWithinBlock.x = t_max.z * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		lSiteLocationWithinBlock.y = t_max.z * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		lSiteLocationWithinBlock.z = t_max.z * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       lSiteLocationWithinBlock,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.z,
			       bCurrentRay);
	      }

	      t_max.z = xyz_Is_1.z
		? t_max.z + t_delta.z
		: t_max.z - t_delta.z;
	    }
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
      

      void ClusterRayTracer::UpdateRayData(const SiteData_t* iSiteData,
				    float ray_t,
				    float ray_segment,
				    Ray* bCurrentRay)
      {
	if (iSiteData->Density < 0.0F)
	{
	  return; // solid voxel
	}

	float palette[3];

	// update the volume rendering of the velocity flow field
	ColPixel::PickColour(iSiteData->Velocity * (float) mDomainStats.velocity_threshold_max_inv,
			     palette);

	UpdateColour(ray_segment, palette, bCurrentRay->VelocityColour);

	if (mVisSettings.mStressType != lb::ShearStress)
	{
	  // update the volume rendering of the von Mises stress flow field
	  float scaled_stress = iSiteData->Stress * (float) mDomainStats.stress_threshold_max_inv;

	  ColPixel::PickColour(scaled_stress, palette);

	  UpdateColour(ray_segment, palette, bCurrentRay->StressColour);
	}

	bCurrentRay->Length += ray_segment;

	if (bCurrentRay->Density >= 0.0F)
	{
	  return;
	}

	bCurrentRay->MinT = ray_t;

	// keep track of the density nearest to the view point
	bCurrentRay->Density = iSiteData->Density;

	// keep track of the stress nearest to the view point
	bCurrentRay->Stress = iSiteData->Stress;		
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
