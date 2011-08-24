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

	lRay.VelocityColour[0] = 0.0F;
	lRay.VelocityColour[1] = 0.0F;
	lRay.VelocityColour[2] = 0.0F;

	lRay.StressColour[0] = 0.0F;
	lRay.StressColour[1] = 0.0F;
	lRay.StressColour[2] = 0.0F;

	lRay.Length = 0.0F;
	lRay.MinT = std::numeric_limits<float>::max();
	lRay.Density = -1.0F;

	Vector3D<bool> lRayInPositiveDirection;
	lRayInPositiveDirection.x = lRay.GetDirection().x > 0.0F;
	lRayInPositiveDirection.y = lRay.GetDirection().y > 0.0F;
	lRayInPositiveDirection.z = lRay.GetDirection().z > 0.0F; 

	TraverseBlocks(&iCluster, lRayInPositiveDirection, lLowerSiteToFirstRayClusterIntersection, &lRay);

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
				     const Vector3D<float>& block_x,
				     const SiteData_t* iSiteData,
				     float t,
				     Ray* bCurrentRay,
				     const Vector3D<bool>& xyz_Is_1)
      {
	site_t i_vec[3];

	i_vec[0] = util::NumericalFunctions::enforceBounds<site_t>((site_t) block_x.x,
								   0,
								   mLatticeData.GetBlockSize() - 1);

	i_vec[1] = util::NumericalFunctions::enforceBounds<site_t>((site_t) block_x.y,
								   0,
								   mLatticeData.GetBlockSize() - 1);

	i_vec[2] = util::NumericalFunctions::enforceBounds<site_t>((site_t) block_x.z,
								   0,
								   mLatticeData.GetBlockSize() - 1);
	

	Vector3D<float> t_max;
	t_max.x = (block_min.x + (float) (xyz_Is_1.x
					  ? i_vec[0] + 1
					  : i_vec[0])) * bCurrentRay->GetInverseDirection().x;

	t_max.y = (block_min.y + (float) (xyz_Is_1.y
					  ? i_vec[1] + 1
					  : i_vec[1])) * bCurrentRay->GetInverseDirection().y;

	t_max.z = (block_min.z + (float) (xyz_Is_1.z
					  ? i_vec[2] + 1
					  : i_vec[2])) * bCurrentRay->GetInverseDirection().z;

	site_t lBlockSizeSquared = mLatticeData.GetBlockSize()*mLatticeData.GetBlockSize();
	site_t lBlockSizeCubed = lBlockSizeSquared*mLatticeData.GetBlockSize();

	site_t i = i_vec[0] * lBlockSizeSquared;
	site_t j = i_vec[1] * mLatticeData.GetBlockSize();
	site_t k = i_vec[2];

	while (true)
	{
	  if (t_max.x < t_max.y)
	  {
	    if (t_max.x < t_max.z)
	    {
	      UpdateRayData(&iSiteData[i+j+k],
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
	      UpdateRayData(&iSiteData[(i + j + k)],
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
	      UpdateRayData(&iSiteData[i + j + k],
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
	      UpdateRayData(&iSiteData[i + j + k],
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

      void ClusterRayTracer::TraverseBlocks(const Cluster* cluster, 
				     const Vector3D<bool>& xyz_Is_1,
				     const Vector3D<float>& ray_dx,
				     Ray *bCurrentRay)
      {

	Vector3D<int> cluster_blocks_vec;
	cluster_blocks_vec.x = cluster->blocksX - 1;
	cluster_blocks_vec.y = cluster->blocksY - 1;
	cluster_blocks_vec.z = cluster->blocksZ - 1;
	int cluster_blocksZ = cluster->blocksZ;
	int cluster_blocksYz = (int) cluster->blocksY * (int) cluster->blocksZ;
	int cluster_blocks = (int) cluster->blocksX * cluster_blocksYz;

	Vector3D<int> i_vec;
	Vector3D<float> block_min;

	float mBlockSizeInverse = 1.0F / mLatticeData.GetBlockSize();
	float mBlockSizeFloat = mLatticeData.GetBlockSize();

	i_vec.x = util::NumericalFunctions::enforceBounds(
	  cluster_blocks_vec.x,
	  0,
	  (int) (mBlockSizeInverse * ray_dx.x));
	
	i_vec.y = util::NumericalFunctions::enforceBounds(
	  cluster_blocks_vec.y,
	  0,
	  (int) (mBlockSizeInverse * ray_dx.y));
	
	i_vec.z = util::NumericalFunctions::enforceBounds(
	  cluster_blocks_vec.z,
	  0,
	  (int) (mBlockSizeInverse * ray_dx.z));

	block_min = Vector3D<float>(i_vec) * mBlockSizeFloat - 
	  Vector3D<float>(ray_dx.x,ray_dx.y,ray_dx.z);
	

	int i = i_vec.x * cluster_blocksYz;
	int j = i_vec.y * cluster_blocksZ;
	int k = i_vec.z;

	Vector3D<float> block_x;
	if (!cluster->SiteData[i + j + k].empty())
	{
	  block_x = block_min * -1.0F;

	  TraverseVoxels(block_min, block_x, &cluster->SiteData[i + j + k][0], 0.0F, bCurrentRay, xyz_Is_1);
	}

	Vector3D <float> t_max;

	t_max.x = (xyz_Is_1.x
		   ? block_min.x + mBlockSizeFloat
		   : block_min.x) * 1.0F * bCurrentRay->GetInverseDirection().x;

	t_max.y = (xyz_Is_1.y
		   ? block_min.y + mBlockSizeFloat
		   : block_min.y) * 1.0F * bCurrentRay->GetInverseDirection().y;

	t_max.z = (xyz_Is_1.z
		   ? block_min.z + mBlockSizeFloat
		   : block_min.z) * 1.0F * bCurrentRay->GetInverseDirection().z;


	Vector3D<float> t_delta = 1.0F * bCurrentRay->GetInverseDirection() * mBlockSizeFloat;
	
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
		block_min.x += mBlockSizeFloat;
	      }
	      else
	      {
		if ( (i -= cluster_blocksYz) < 0)
		  return;
		block_min.x -= mBlockSizeFloat;
	      }

	      if (!cluster->SiteData[i + j + k].empty())
	      {
		block_x.x = t_max.x * bCurrentRay->GetDirection().x - block_min.x;
		block_x.y = t_max.x * bCurrentRay->GetDirection().y - block_min.y;
		block_x.z = t_max.x * bCurrentRay->GetDirection().z - block_min.z;

		TraverseVoxels(block_min,
			       block_x,
			       &cluster->SiteData[i + j + k][0],
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
		block_min.z += mBlockSizeFloat;
	      }
	      else
	      {
		if (--k < 0)
		  return;
		block_min.z -= mBlockSizeFloat;
	      }

	      if (!cluster->SiteData[i + j + k].empty())
	      {
		block_x.x = t_max.z * bCurrentRay->GetDirection().x - block_min.x;
		block_x.y = t_max.z * bCurrentRay->GetDirection().y - block_min.y;
		block_x.z = t_max.z * bCurrentRay->GetDirection().z - block_min.z;

		TraverseVoxels(block_min,
			       block_x,
			       &cluster->SiteData[i + j + k][0],
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
		block_min.y += mBlockSizeFloat;
	      }
	      else
	      {
		if ( (j -= cluster_blocksZ) < 0)
		  return;
		block_min.y -= mBlockSizeFloat;
	      }

	      if (!cluster->SiteData[i + j + k].empty())
	      {
		block_x.x = t_max.y * bCurrentRay->GetDirection().x - block_min.x;
		block_x.y = t_max.y * bCurrentRay->GetDirection().y - block_min.y;
		block_x.z = t_max.y * bCurrentRay->GetDirection().z - block_min.z;

		TraverseVoxels(block_min,
			       block_x,
			       &cluster->SiteData[i + j + k][0],
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
		block_min.z += mBlockSizeFloat;
	      }
	      else
	      {
		if (--k < 0)
		  return;
		block_min.z -= mBlockSizeFloat;
	      }

	      if (!cluster->SiteData[i + j + k].empty())
	      {
		block_x.x = t_max.z * bCurrentRay->GetDirection().x - block_min.x;
		block_x.y = t_max.z * bCurrentRay->GetDirection().y - block_min.y;
		block_x.z = t_max.z * bCurrentRay->GetDirection().z - block_min.z;

		TraverseVoxels(block_min,
			       block_x,
			       &cluster->SiteData[i + j + k][0],
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
