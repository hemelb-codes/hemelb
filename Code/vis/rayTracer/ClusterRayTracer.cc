//#define NDEBUG;
#include <cmath> 
#include <cassert>
#include <iostream>

#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/ClusterRayTracer.h"
#include "vis/rayTracer/Ray.h"

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

      void ClusterRayTracer::TraverseVoxels(const Vector3D<float>& block_min,
					    const Vector3D<float>& iLocationInBlock,
					    const SiteData_t* iStartOfSiteDataForBlock,
					    float t,
					    Ray* bCurrentRay,
					    const Vector3D<bool>& xyz_Is_1)
      {
	Vector3D<site_t> lLocationInBlock = EnforceBlockBounds(iLocationInBlock);

	Vector3D<float> t_max;
	t_max.x = (block_min.x + (float) (xyz_Is_1.x
					  ? lLocationInBlock.x + 1
					  : lLocationInBlock.x)) * bCurrentRay->GetInverseDirection().x;

	t_max.y = (block_min.y + (float) (xyz_Is_1.y
					  ? lLocationInBlock.y + 1
					  : lLocationInBlock.y)) * bCurrentRay->GetInverseDirection().y;

	t_max.z = (block_min.z + (float) (xyz_Is_1.z
					  ? lLocationInBlock.z + 1
					  : lLocationInBlock.z)) * bCurrentRay->GetInverseDirection().z;

	site_t lBlockSizeSquared = mLatticeData.GetBlockSize()*mLatticeData.GetBlockSize();
	site_t lBlockSizeCubed = lBlockSizeSquared*mLatticeData.GetBlockSize();

	site_t i = lLocationInBlock.x * lBlockSizeSquared;
	site_t j = lLocationInBlock.y * mLatticeData.GetBlockSize();
	site_t k = lLocationInBlock.z;

	while (true)
	{
	  if (t_max.x < t_max.y)
	  {
	    if (t_max.x < t_max.z)
	    {
	      UpdateRayData(&iStartOfSiteDataForBlock[i+j+k],
			    t,
			    t_max.x - t,
			    bCurrentRay);

	      if (xyz_Is_1.x)
	      {
		if ( (i += lBlockSizeSquared) >= lBlockSizeCubed)
		{
		  return;
		}
		t = t_max.x;
		t_max.x += 1.0F * bCurrentRay->GetInverseDirection().x;
	      }
	      else
	      {
		if (i < lBlockSizeSquared)
		{
		  return;
		}
		else
		{
		  i -= lBlockSizeSquared;
		}
		t = t_max.x;
		t_max.x -= 1.0F * bCurrentRay->GetInverseDirection().x;
	      }
	    }
	    else
	    {
	      UpdateRayData(&iStartOfSiteDataForBlock[(i + j + k)],
			    t,
			    t_max.z - t,
			    bCurrentRay);

	      if (xyz_Is_1.z)
	      {
		if (++k >= mLatticeData.GetBlockSize())
		{
		  return;
		}
		t = t_max.z;
		t_max.z += 1.0F * bCurrentRay->GetInverseDirection().z;
	      }
	      else
	      {
		if (k == 0)
		{
		  return;
		}
		else
		{
		  --k;
		}
		t = t_max.z;
		t_max.z -= 1.0F * bCurrentRay->GetInverseDirection().z;
	      }
	    }
	  }
	  else
	  {
	    if (t_max.y < t_max.z)
	    {
	      UpdateRayData(&iStartOfSiteDataForBlock[i + j + k],
			    t,
			    t_max.y - t,
			    bCurrentRay);

	      if (xyz_Is_1.y)
	      {
		if ( (j += mLatticeData.GetBlockSize()) >= lBlockSizeSquared)
		{
		  return;
		}
		t = t_max.y;
		t_max.y += 1.0F * bCurrentRay->GetInverseDirection().y;
	      }
	      else
	      {
		if (j < mLatticeData.GetBlockSize())
		{
		  return;
		}
		else
		{
		  j -= mLatticeData.GetBlockSize();
		}
		t = t_max.y;
		t_max.y -= 1.0F * bCurrentRay->GetInverseDirection().y;
	      }
	    }
	    else
	    {
	      UpdateRayData(&iStartOfSiteDataForBlock[i + j + k],
			    t,
			    t_max.z - t,
			    bCurrentRay);

	      if (xyz_Is_1.z)
	      {
		if (++k >= mLatticeData.GetBlockSize())
		{
		  return;
		}
		t = t_max.z;
		t_max.z += 1.0F * bCurrentRay->GetInverseDirection().z;
	      }
	      else
	      {
		if (k == 0)
		{
		  return;
		}
		else
		{
		  --k;
		}
		t = t_max.z;
		t_max.z -= 1.0F * bCurrentRay->GetInverseDirection().z;
	      }
	    }
	  }
	}
      }

      Vector3D<site_t> ClusterRayTracer::
      EnforceBlockBounds(const Vector3D<float>& iUnboundLocation)
      {
	Vector3D<site_t> lLocationInBlock;
	lLocationInBlock.x = util::NumericalFunctions::
	  enforceBounds<site_t>((site_t) iUnboundLocation.x,
				0,
				mLatticeData.GetBlockSize() - 1);
	  
	lLocationInBlock.y  = util::NumericalFunctions::
	  enforceBounds<site_t>((site_t) iUnboundLocation.y,
				0,
				mLatticeData.GetBlockSize() - 1);

	lLocationInBlock.z = util::NumericalFunctions::
	  enforceBounds<site_t>((site_t) iUnboundLocation.z,
				0,
				mLatticeData.GetBlockSize() - 1);

	return lLocationInBlock;
      }
 

      void ClusterRayTracer::TraverseBlocks(const Cluster& iCluster, 
				     const Vector3D<bool>& xyz_Is_1,
				     const Vector3D<float>& iLowerSiteToFirstRayClusterIntersection,
				     Ray *bCurrentRay)
      {
	int cluster_blocksZ = iCluster.blocksZ;
	int cluster_blocksYz = (int) iCluster.blocksY * (int) iCluster.blocksZ;
	int cluster_blocks = (int) iCluster.blocksX * cluster_blocksYz;

	Vector3D<float> block_min;

	float lBlockSizeFloat = mLatticeData.GetBlockSize();

 	Vector3D<unsigned int> lBlockCoordinatesOfFirstIntersectionBlock = 
	  GetBlockCoordinatesOfFirstIntersectionBlock(iCluster, iLowerSiteToFirstRayClusterIntersection);
	
	Vector3D<float> lFirstIntersectionToBlockLowerSite = 
	  Vector3D<float>(lBlockCoordinatesOfFirstIntersectionBlock) * lBlockSizeFloat - 
	  iLowerSiteToFirstRayClusterIntersection;
	
	int i = lBlockCoordinatesOfFirstIntersectionBlock.x * cluster_blocksYz;
	int j = lBlockCoordinatesOfFirstIntersectionBlock.y * cluster_blocksZ;
	int k = lBlockCoordinatesOfFirstIntersectionBlock.z; 

	//unsigned int lBlockId = iCluster.GetBlockIdFrom3DBlockLocation(lBlockCoordinatesOfFirstIntersectionBlock);

	Vector3D<float> block_x;
	if (!iCluster.SiteData[i + j + k].empty())
	{
	  block_x = lFirstIntersectionToBlockLowerSite * -1.0F;

	  TraverseVoxels(lFirstIntersectionToBlockLowerSite, block_x, &iCluster.SiteData[i + j + k][0], 0.0F, bCurrentRay, xyz_Is_1);
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
		block_x.x = t_max.x * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		block_x.y = t_max.x * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		block_x.z = t_max.x * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       block_x,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.x,
			       bCurrentRay,
			       xyz_Is_1);
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
		block_x.x = t_max.z * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		block_x.y = t_max.z * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		block_x.z = t_max.z * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       block_x,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.z,
			       bCurrentRay,
			       xyz_Is_1);
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
		block_x.x = t_max.y * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		block_x.y = t_max.y * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		block_x.z = t_max.y * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       block_x,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.y,
			       bCurrentRay,
			       xyz_Is_1);
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
		block_x.x = t_max.z * bCurrentRay->GetDirection().x - lFirstIntersectionToBlockLowerSite.x;
		block_x.y = t_max.z * bCurrentRay->GetDirection().y - lFirstIntersectionToBlockLowerSite.y;
		block_x.z = t_max.z * bCurrentRay->GetDirection().z - lFirstIntersectionToBlockLowerSite.z;

		TraverseVoxels(lFirstIntersectionToBlockLowerSite,
			       block_x,
			       &iCluster.SiteData[i + j + k][0],
			       t_max.z,
			       bCurrentRay,
			       xyz_Is_1);
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
	Vector3D<unsigned int> lBlockCoordinatesOFFirstIntersectionBlock;

	//Perform the truncated division and ensure that the 
	//coordinates are valid to allow for numerical errors
	
	lBlockCoordinatesOFFirstIntersectionBlock.x = (unsigned int)
	  util::NumericalFunctions::enforceBounds(
	  iCluster.blocksX - 1,
	  0,
	  (int) (1.0F / mLatticeData.GetBlockSize() * iLowerSiteToFirstRayClusterIntersection.x));
	
	lBlockCoordinatesOFFirstIntersectionBlock.y = (unsigned int)
	  util::NumericalFunctions::enforceBounds(
	  iCluster.blocksY - 1,
	  0,
	  (int) (1.0F / mLatticeData.GetBlockSize() * iLowerSiteToFirstRayClusterIntersection.y));
	
	lBlockCoordinatesOFFirstIntersectionBlock.z = (unsigned int)
	  util::NumericalFunctions::enforceBounds(
	  iCluster.blocksZ - 1,
	  0,
	  (int) (1.0F / mLatticeData.GetBlockSize() * iLowerSiteToFirstRayClusterIntersection.z));

	return lBlockCoordinatesOFFirstIntersectionBlock;
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
