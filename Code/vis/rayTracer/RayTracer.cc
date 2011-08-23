//#define NDEBUG;
#include <assert.h>


#include <math.h>
#include <stdlib.h>
#include <vector>
#include <limits>

#include "debug/Debugger.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h" 
#include "vis/Vector3D.h"
#include "vis/XYCoordinates.h"
#include "vis/rayTracer/RayTracer.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      // TODO RENAME THIS FUNCTION
      void RayTracer::AABBvsRay(const AABB* aabb,
				const Vector3D<float>& inverseDirection,
				const Vector3D<bool>& xyzComponentIsPositive,
				float* t_near,
				float* t_far)
      {
	float tx0 = (xyzComponentIsPositive.x
		     ? aabb->acc_2
		     : aabb->acc_1) * inverseDirection.x;
	float tx1 = (xyzComponentIsPositive.x
		     ? aabb->acc_1
		     : aabb->acc_2) * inverseDirection.x;

	float ty0 = (xyzComponentIsPositive.y
		     ? aabb->acc_4
		     : aabb->acc_3) * inverseDirection.y;
	float ty1 = (xyzComponentIsPositive.y
		     ? aabb->acc_3
		     : aabb->acc_4) * inverseDirection.y;

	float tz0 = (xyzComponentIsPositive.z
		     ? aabb->acc_6
		     : aabb->acc_5) * inverseDirection.z;
	float tz1 = (xyzComponentIsPositive.z
		     ? aabb->acc_5
		     : aabb->acc_6) * inverseDirection.z;

	*t_near = fmaxf(tx0, fmaxf(ty0, tz0));
	*t_far = fminf(tx1, fminf(ty1, tz1));
      }

      void RayTracer::BuildClusters()
      {
	mClusterBuilder = new ClusterBuilder(mLatDat);
	mClusterBuilder->BuildClusters();
      }

      /**
       * Update a colour vector to include a section of known length through a
       * solid of known colour.
       *
       * @param dt
       * @param palette
       * @param col
       */
      void RayTracer::UpdateColour(float dt, const float palette[3], float col[3])
      {
	for (int ii = 0; ii < 3; ++ii)
	{
	  col[ii] += dt * palette[ii];
	}
      }

      void RayTracer::UpdateRayData(const SiteData_t* iSiteData,
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
	ColPixel::PickColour(iSiteData->Velocity * (float) mDomainStats->velocity_threshold_max_inv,
			     palette);

	UpdateColour(ray_segment, palette, bCurrentRay->VelocityColour);

	if (mVisSettings->mStressType != lb::ShearStress)
	{
	  // update the volume rendering of the von Mises stress flow field
	  float scaled_stress = iSiteData->Stress * (float) mDomainStats->stress_threshold_max_inv;

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

      void RayTracer::TraverseVoxels(const Vector3D<float>& block_min,
				     const Vector3D<float>& block_x,
				     const SiteData_t* iSiteData,
				     float t,
				     Ray* bCurrentRay,
				     const Vector3D<bool>& xyz_Is_1)
      {
	site_t i_vec[3];

	i_vec[0] = util::NumericalFunctions::enforceBounds<site_t>((site_t) block_x.x,
								     0,
								     block_size_1);

	i_vec[1] = util::NumericalFunctions::enforceBounds<site_t>((site_t) block_x.y,
								     0,
								     block_size_1);

	i_vec[2] = util::NumericalFunctions::enforceBounds<site_t>((site_t) block_x.z,
								     0,
								     block_size_1);
	

	Vector3D<float> t_max;
	t_max.x = (block_min.x + (float) (xyz_Is_1.x
					      ? i_vec[0] + 1
					      : i_vec[0])) * bCurrentRay->InverseDirection.x;

	t_max.y = (block_min.y + (float) (xyz_Is_1.y
					      ? i_vec[1] + 1
					      : i_vec[1])) * bCurrentRay->InverseDirection.y;

	t_max.z = (block_min.z + (float) (xyz_Is_1.z
					      ? i_vec[2] + 1
					      : i_vec[2])) * bCurrentRay->InverseDirection.z;


	site_t i = i_vec[0] * block_size2;
	site_t j = i_vec[1] * mLatDat->GetBlockSize();
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
		if ( (i += block_size2) >= block_size3)
		{
		  return;
		}
		t = t_max.x;
		t_max.x += bCurrentRay->InverseDirection.x;
	      }
	      else
	      {
		if (i < block_size2)
		{
		  return;
		}
		else
		{
		  i -= block_size2;
		}
		t = t_max.x;
		t_max.x -= bCurrentRay->InverseDirection.x;
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
		if (++k >= mLatDat->GetBlockSize())
		{
		  return;
		}
		t = t_max.z;
		t_max.z += bCurrentRay->InverseDirection.z;
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
		t_max.z -= bCurrentRay->InverseDirection.z;
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
		if ( (j += mLatDat->GetBlockSize()) >= block_size2)
		{
		  return;
		}
		t = t_max.y;
		t_max.y += bCurrentRay->InverseDirection.y;
	      }
	      else
	      {
		if (j < mLatDat->GetBlockSize())
		{
		  return;
		}
		else
		{
		  j -= mLatDat->GetBlockSize();
		}
		t = t_max.y;
		t_max.y -= bCurrentRay->InverseDirection.y;
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
		if (++k >= mLatDat->GetBlockSize())
		{
		  return;
		}
		t = t_max.z;
		t_max.z += bCurrentRay->InverseDirection.z;
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
		t_max.z -= bCurrentRay->InverseDirection.z;
	      }
	    }
	  }
	}
      }

      void RayTracer::TraverseBlocks(const Cluster* cluster, 
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
		   : block_min.x) * bCurrentRay->InverseDirection.x;

	t_max.y = (xyz_Is_1.y
		    ? block_min.y + mBlockSizeFloat
		   : block_min.y) * bCurrentRay->InverseDirection.y;

	t_max.z = (xyz_Is_1.z
		   ? block_min.z + mBlockSizeFloat
		   : block_min.z) * bCurrentRay->InverseDirection.z;


	Vector3D<float> t_delta = bCurrentRay->InverseDirection * mBlockSizeFloat;
	
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
		block_x.x = t_max.x * bCurrentRay->Direction.x - block_min.x;
		block_x.y = t_max.x * bCurrentRay->Direction.y - block_min.y;
		block_x.z = t_max.x * bCurrentRay->Direction.z - block_min.z;

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
		block_x.x = t_max.z * bCurrentRay->Direction.x - block_min.x;
		block_x.y = t_max.z * bCurrentRay->Direction.y - block_min.y;
		block_x.z = t_max.z * bCurrentRay->Direction.z - block_min.z;

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
		block_x.x = t_max.y * bCurrentRay->Direction.x - block_min.x;
		block_x.y = t_max.y * bCurrentRay->Direction.y - block_min.y;
		block_x.z = t_max.y * bCurrentRay->Direction.z - block_min.z;

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
		block_x.x = t_max.z * bCurrentRay->Direction.x - block_min.x;
		block_x.y = t_max.z * bCurrentRay->Direction.y - block_min.y;
		block_x.z = t_max.z * bCurrentRay->Direction.z - block_min.z;

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

      RayTracer::RayTracer(const geometry::LatticeData* iLatDat,
			   const DomainStats* iDomainStats,
			   Screen* iScreen,
			   Viewpoint* iViewpoint,
			   VisSettings* iVisSettings) :
  	mLatDat(iLatDat), mDomainStats(iDomainStats), mScreen(iScreen), mViewpoint(iViewpoint),
	mVisSettings(iVisSettings), mBlockSizeFloat((float) mLatDat->GetBlockSize()),
	mBlockSizeInverse(1.F / mBlockSizeFloat), block_size2(mLatDat->GetBlockSize()
							      * mLatDat->GetBlockSize()), block_size3(mLatDat->GetBlockSize() * block_size2),
	block_size_1(mLatDat->GetBlockSize() - 1), blocksYz(mLatDat->GetYBlockCount()
							    * mLatDat->GetZBlockCount())
      {
	mClusterBuilder = NULL;
      }

      void RayTracer::Render()
      {
	const Vector3D<float>& projectedUnitX = mScreen->GetUnitVectorProjectionX();
	const Vector3D<float>& projectedUnitY = mScreen->GetUnitVectorProjectionY();

	for (unsigned int clusterId = 0; clusterId < mClusterBuilder->GetClusters().size(); clusterId++)
	{
	  RenderCluster(*mClusterBuilder->GetClusters()[clusterId]);
	}
      }

      void RayTracer::RenderCluster(const Cluster& iCluster)
      {
	Vector3D <float> lLowerSiteCordinatesOfClusterRelativeToViewpoint = 
	  iCluster.minBlock - mViewpoint->GetViewpointCentreLocation();
	  
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
	
	XYCoordinates<int> lSubImageLowerLeftPixelCoordinates = 
	  mScreen->TransformScreenToPixelCoordinates<int> (lSubImageLowerLeft);

	XYCoordinates<int> lSubImageUpperRightPixelCoordinates = 
	  mScreen->TransformScreenToPixelCoordinates<int> (lSubImageUpperRight);

	// If the entire sub-image is off the screen, continue to the next cluster.
	if (lSubImageLowerLeftPixelCoordinates.x >= mScreen->GetPixelsX() || lSubImageUpperRightPixelCoordinates.x < 0 || lSubImageLowerLeftPixelCoordinates.y
	    >= mScreen->GetPixelsY() || lSubImageUpperRightPixelCoordinates.y < 0)
	{
	  return;
	}

	// Crop the sub-image to the screen.
	lSubImageLowerLeftPixelCoordinates.x = util::NumericalFunctions::
	  max(lSubImageLowerLeftPixelCoordinates.x, 0);
	
	lSubImageUpperRightPixelCoordinates.x = util::NumericalFunctions::
	  min(lSubImageUpperRightPixelCoordinates.x, mScreen->GetPixelsX() - 1);
	
	lSubImageLowerLeftPixelCoordinates.y = util::NumericalFunctions::
	  max(lSubImageLowerLeftPixelCoordinates.y, 0);

	lSubImageUpperRightPixelCoordinates.y = util::NumericalFunctions::
	  min(lSubImageUpperRightPixelCoordinates.y, mScreen->GetPixelsY() - 1);

	AABB aabb;
	aabb.acc_1 = iCluster.maxSite.x - mViewpoint->GetViewpointCentreLocation().x;
	aabb.acc_2 = iCluster.minSite.x - mViewpoint->GetViewpointCentreLocation().x;
	aabb.acc_3 = iCluster.maxSite.y - mViewpoint->GetViewpointCentreLocation().y;
	aabb.acc_4 = iCluster.minSite.y - mViewpoint->GetViewpointCentreLocation().y;
	aabb.acc_5 = iCluster.maxSite.z - mViewpoint->GetViewpointCentreLocation().z;
	aabb.acc_6 = iCluster.minSite.z - mViewpoint->GetViewpointCentreLocation().z;

	Vector3D<float> par3;
	const Vector3D<float>& vtx = mScreen->GetVtx();
	  
	par3 = vtx + mScreen->GetUnitVectorProjectionX() * (float) lSubImageLowerLeftPixelCoordinates.x
	  + mScreen->GetUnitVectorProjectionY() * (float) lSubImageLowerLeftPixelCoordinates.y;
	  

	  for (int subImageX = lSubImageLowerLeftPixelCoordinates.x; subImageX <= lSubImageUpperRightPixelCoordinates.x; ++subImageX)
	  {
	    Vector3D<float> lRayDirection = par3;

	    for (int subImageY = lSubImageLowerLeftPixelCoordinates.y; subImageY <= lSubImageUpperRightPixelCoordinates.y; ++subImageY)
	    {
	      Ray lRay;

	      lRay.Direction.x = lRayDirection.x;
	      lRay.Direction.y = lRayDirection.y;
	      lRay.Direction.z = lRayDirection.z;

	      float lInverseDirectionMagnitude = 1.0F / sqrtf(lRayDirection.x * lRayDirection.x
							      + lRayDirection.y * lRayDirection.y + lRayDirection.z * lRayDirection.z);

	      lRay.Direction.x *= lInverseDirectionMagnitude;
	      lRay.Direction.y *= lInverseDirectionMagnitude;
	      lRay.Direction.z *= lInverseDirectionMagnitude;
	       
	      lRay.InverseDirection.x = 1.0F / lRay.Direction.x;
	      lRay.InverseDirection.y = 1.0F / lRay.Direction.y;
	      lRay.InverseDirection.z = 1.0F / lRay.Direction.z;

	      Vector3D<bool> lRayInPositiveDirection;
	      lRayInPositiveDirection.x = lRay.Direction.x > 0.0F;
	      lRayInPositiveDirection.y = lRay.Direction.y > 0.0F;
	      lRayInPositiveDirection.z = lRay.Direction.z > 0.0F;

	      lRayDirection.x += mScreen->GetUnitVectorProjectionY().x;
	      lRayDirection.y += mScreen->GetUnitVectorProjectionY().y;
	      lRayDirection.z += mScreen->GetUnitVectorProjectionY().z;

	      float t_near, t_far;
	      AABBvsRay(&aabb, lRay.InverseDirection, lRayInPositiveDirection, &t_near, &t_far);

	      if (t_near > t_far)
	      {
		continue;
	      }

	      Vector3D <float> ray_dx = t_near * lRay.Direction -
		lLowerSiteCordinatesOfClusterRelativeToViewpoint;

	      lRay.VelocityColour[0] = 0.0F;
	      lRay.VelocityColour[1] = 0.0F;
	      lRay.VelocityColour[2] = 0.0F;

	      lRay.StressColour[0] = 0.0F;
	      lRay.StressColour[1] = 0.0F;
	      lRay.StressColour[2] = 0.0F;

	      lRay.Length = 0.0F;
	      lRay.MinT = std::numeric_limits<float>::max();
	      lRay.Density = -1.0F;

	      TraverseBlocks(&iCluster, lRayInPositiveDirection, ray_dx, &lRay);

	      if (lRay.MinT == std::numeric_limits<float>::max())
	      {
		continue;
	      }

	      ColPixel col_pixel(subImageX, subImageY, lRay.MinT + t_near, lRay.Length, (lRay.Density
											 - (float) mDomainStats->density_threshold_min)
				 * (float) mDomainStats->density_threshold_minmax_inv, lRay.Stress
				 != std::numeric_limits<float>::max()
				 ? lRay.Stress * (float) mDomainStats->stress_threshold_max_inv
				 : std::numeric_limits<float>::max(), lRay.VelocityColour, lRay.StressColour);

	      mScreen->AddPixel(&col_pixel, mVisSettings);
	    }
	    par3+=mScreen->GetUnitVectorProjectionX();
	  }
      }

      void RayTracer::UpdateSubImageExtentForCorner
      (const Vector3D<float>& iCorner,
       XYCoordinates<float>& ioSubImageLowerLeft,
       XYCoordinates<float>& ioSubImageUpperRight)
      {
	XYCoordinates<float> lCornerProjection = mViewpoint->FlatProject(iCorner);

	XYCoordinates<float>::UpdateMinXYCoordinates(ioSubImageLowerLeft, lCornerProjection);
	XYCoordinates<float>::UpdateMaxXYCoordinates(ioSubImageUpperRight, lCornerProjection);
      }

      void RayTracer::UpdateClusterVoxel(site_t i,
					 distribn_t density,
					 distribn_t velocity,
					 distribn_t stress)
      {
	assert(static_cast<site_t>(static_cast<unsigned int>(i)) == i);
	
	
	mClusterBuilder->GetClusterVoxelDataPointer(i)->Density =
	  (float) density;
	mClusterBuilder->GetClusterVoxelDataPointer(i)->Velocity =
	  (float) velocity;
	mClusterBuilder->GetClusterVoxelDataPointer(i)->Stress =
	  (float) stress;
      }

      RayTracer::~RayTracer()
      {
	delete mClusterBuilder;
      }
    }
  }
}
