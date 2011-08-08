#include <math.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <assert.h>

#include "debug/Debugger.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "vis/rayTracer/Location.h"
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
				const float inverseDirection[3],
				const bool xyzComponentIsPositive[3],
				float* t_near,
				float* t_far)
      {
	float tx0 = (xyzComponentIsPositive[0]
		     ? aabb->acc_2
		     : aabb->acc_1) * inverseDirection[0];
	float tx1 = (xyzComponentIsPositive[0]
		     ? aabb->acc_1
		     : aabb->acc_2) * inverseDirection[0];
	float ty0 = (xyzComponentIsPositive[1]
		     ? aabb->acc_4
		     : aabb->acc_3) * inverseDirection[1];
	float ty1 = (xyzComponentIsPositive[1]
		     ? aabb->acc_3
		     : aabb->acc_4) * inverseDirection[1];
	float tz0 = (xyzComponentIsPositive[2]
		     ? aabb->acc_6
		     : aabb->acc_5) * inverseDirection[2];
	float tz1 = (xyzComponentIsPositive[2]
		     ? aabb->acc_5
		     : aabb->acc_6) * inverseDirection[2];

	*t_near = fmaxf(tx0, fmaxf(ty0, tz0));
	*t_far = fminf(tx1, fminf(ty1, tz1));
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

      void RayTracer::UpdateRayData(const float flow_field[3],
				    float ray_t,
				    float ray_segment,
				    Ray* bCurrentRay)
      {
	if (flow_field[0] < 0.0F)
	{
	  return; // solid voxel
	}

	float palette[3];

	// update the volume rendering of the velocity flow field
	ColPixel::PickColour(flow_field[1] * (float) mDomainStats->velocity_threshold_max_inv,
			     palette);

	UpdateColour(ray_segment, palette, bCurrentRay->VelocityColour);

	if (mVisSettings->mStressType != lb::ShearStress)
	{
	  // update the volume rendering of the von Mises stress flow field
	  float scaled_stress = flow_field[2] * (float) mDomainStats->stress_threshold_max_inv;

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
	bCurrentRay->Density = flow_field[0];

	// keep track of the stress nearest to the view point
	bCurrentRay->Stress = flow_field[2];		
      }

      void RayTracer::TraverseVoxels(const float block_min[3],
				     const float block_x[3],
				     const float voxel_flow_field[],
				     float t,
				     Ray* bCurrentRay,
				     const bool xyz_is_1[3])
      {
	site_t i_vec[3];

	for (int i = 0; i < 3; i++)
	{
	  i_vec[i] = util::NumericalFunctions::enforceBounds<site_t>((site_t) block_x[i],
								     0,
								     block_size_1);
	}

	float t_max[3];
	for (int i = 0; i < 3; i++)
	{
	  t_max[i] = (block_min[i] + (float) (xyz_is_1[i]
					      ? i_vec[i] + 1
					      : i_vec[i])) * bCurrentRay->InverseDirection[i];
	}

	site_t i = i_vec[0] * block_size2;
	site_t j = i_vec[1] * mLatDat->GetBlockSize();
	site_t k = i_vec[2];

	while (true)
	{
	  if (t_max[0] < t_max[1])
	  {
	    if (t_max[0] < t_max[2])
	    {
	      UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
			    t,
			    t_max[0] - t,
			    bCurrentRay);

	      if (xyz_is_1[0])
	      {
		if ( (i += block_size2) >= block_size3)
		{
		  return;
		}
		t = t_max[0];
		t_max[0] += bCurrentRay->InverseDirection[0];
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
		t = t_max[0];
		t_max[0] -= bCurrentRay->InverseDirection[0];
	      }
	    }
	    else
	    {
	      UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
			    t,
			    t_max[2] - t,
			    bCurrentRay);

	      if (xyz_is_1[2])
	      {
		if (++k >= mLatDat->GetBlockSize())
		{
		  return;
		}
		t = t_max[2];
		t_max[2] += bCurrentRay->InverseDirection[2];
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
		t = t_max[2];
		t_max[2] -= bCurrentRay->InverseDirection[2];
	      }
	    }
	  }
	  else
	  {
	    if (t_max[1] < t_max[2])
	    {
	      UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
			    t,
			    t_max[1] - t,
			    bCurrentRay);

	      if (xyz_is_1[1])
	      {
		if ( (j += mLatDat->GetBlockSize()) >= block_size2)
		{
		  return;
		}
		t = t_max[1];
		t_max[1] += bCurrentRay->InverseDirection[1];
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
		t = t_max[1];
		t_max[1] -= bCurrentRay->InverseDirection[1];
	      }
	    }
	    else
	    {
	      UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
			    t,
			    t_max[2] - t,
			    bCurrentRay);

	      if (xyz_is_1[2])
	      {
		if (++k >= mLatDat->GetBlockSize())
		{
		  return;
		}
		t = t_max[2];
		t_max[2] += bCurrentRay->InverseDirection[2];
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
		t = t_max[2];
		t_max[2] -= bCurrentRay->InverseDirection[2];
	      }
	    }
	  }
	}
      }

      void RayTracer::TraverseBlocks(const Cluster* cluster,
				     const bool xyz_Is_1[3],
				     const float ray_dx[3],
				     float **block_flow_field,
				     Ray *bCurrentRay)
      {

	int cluster_blocks_vec[3];
	cluster_blocks_vec[0] = cluster->blocksX - 1;
	cluster_blocks_vec[1] = cluster->blocksY - 1;
	cluster_blocks_vec[2] = cluster->blocksZ - 1;
	int cluster_blocksZ = cluster->blocksZ;
	int cluster_blocksYz = (int) cluster->blocksY * (int) cluster->blocksZ;
	int cluster_blocks = (int) cluster->blocksX * cluster_blocksYz;

	int i_vec[3];
	float block_min[3];

	for (int l = 0; l < 3; l++)
	{
	  i_vec[l] = util::NumericalFunctions::enforceBounds(cluster_blocks_vec[l],
							     0,
							     (int) (mBlockSizeInverse * ray_dx[l]));
	  block_min[l] = (float) i_vec[l] * mBlockSizeFloat - ray_dx[l];
	}

	int i = i_vec[0] * cluster_blocksYz;
	int j = i_vec[1] * cluster_blocksZ;
	int k = i_vec[2];

	float block_x[3];
	if (block_flow_field[i + j + k] != NULL)
	{
	  block_x[0] = -block_min[0];
	  block_x[1] = -block_min[1];
	  block_x[2] = -block_min[2];

	  TraverseVoxels(block_min, block_x, block_flow_field[i + j + k], 0.0F, bCurrentRay, xyz_Is_1);
	}

	float t_max[3];
	float t_delta[3];
	for (int l = 0; l < 3; l++)
	{
	  t_max[l] = (xyz_Is_1[l]
		      ? block_min[l] + mBlockSizeFloat
		      : block_min[l]) * bCurrentRay->InverseDirection[l];
	  t_delta[l] = mBlockSizeFloat * bCurrentRay->InverseDirection[l];
	}

	while (true)
	{
	  if (t_max[0] < t_max[1])
	  {
	    if (t_max[0] < t_max[2])
	    {
	      if (xyz_Is_1[0])
	      {
		if ( (i += cluster_blocksYz) >= cluster_blocks)
		  return;
		block_min[0] += mBlockSizeFloat;
	      }
	      else
	      {
		if ( (i -= cluster_blocksYz) < 0)
		  return;
		block_min[0] -= mBlockSizeFloat;
	      }

	      if (block_flow_field[i + j + k] != NULL)
	      {
		block_x[0] = t_max[0] * bCurrentRay->Direction[0] - block_min[0];
		block_x[1] = t_max[0] * bCurrentRay->Direction[1] - block_min[1];
		block_x[2] = t_max[0] * bCurrentRay->Direction[2] - block_min[2];

		TraverseVoxels(block_min,
			       block_x,
			       block_flow_field[i + j + k],
			       t_max[0],
			       bCurrentRay,
			       xyz_Is_1);
	      }

	      t_max[0] = xyz_Is_1[0]
		? t_max[0] + t_delta[0]
		: t_max[0] - t_delta[0];
	    }
	    else
	    {
	      if (xyz_Is_1[2])
	      {
		if (++k >= cluster_blocksZ)
		  return;
		block_min[2] += mBlockSizeFloat;
	      }
	      else
	      {
		if (--k < 0)
		  return;
		block_min[2] -= mBlockSizeFloat;
	      }

	      if (block_flow_field[i + j + k] != NULL)
	      {
		block_x[0] = t_max[2] * bCurrentRay->Direction[0] - block_min[0];
		block_x[1] = t_max[2] * bCurrentRay->Direction[1] - block_min[1];
		block_x[2] = t_max[2] * bCurrentRay->Direction[2] - block_min[2];

		TraverseVoxels(block_min,
			       block_x,
			       block_flow_field[i + j + k],
			       t_max[2],
			       bCurrentRay,
			       xyz_Is_1);
	      }

	      t_max[2] = xyz_Is_1[2]
		? t_max[2] + t_delta[2]
		: t_max[2] - t_delta[2];
	    }
	  }
	  else
	  {
	    if (t_max[1] < t_max[2])
	    {
	      if (xyz_Is_1[1])
	      {
		if ( (j += cluster_blocksZ) >= cluster_blocksYz)
		  return;
		block_min[1] += mBlockSizeFloat;
	      }
	      else
	      {
		if ( (j -= cluster_blocksZ) < 0)
		  return;
		block_min[1] -= mBlockSizeFloat;
	      }

	      if (block_flow_field[i + j + k] != NULL)
	      {
		block_x[0] = t_max[1] * bCurrentRay->Direction[0] - block_min[0];
		block_x[1] = t_max[1] * bCurrentRay->Direction[1] - block_min[1];
		block_x[2] = t_max[1] * bCurrentRay->Direction[2] - block_min[2];

		TraverseVoxels(block_min,
			       block_x,
			       block_flow_field[i + j + k],
			       t_max[1],
			       bCurrentRay,
			       xyz_Is_1);
	      }

	      t_max[1] = xyz_Is_1[1]
		? t_max[1] + t_delta[1]
		: t_max[1] - t_delta[1];
	    }
	    else
	    {
	      if (xyz_Is_1[2])
	      {
		if (++k >= cluster_blocksZ)
		  return;
		block_min[2] += mBlockSizeFloat;
	      }
	      else
	      {
		if (--k < 0)
		  return;
		block_min[2] -= mBlockSizeFloat;
	      }

	      if (block_flow_field[i + j + k] != NULL)
	      {
		block_x[0] = t_max[2] * bCurrentRay->Direction[0] - block_min[0];
		block_x[1] = t_max[2] * bCurrentRay->Direction[1] - block_min[1];
		block_x[2] = t_max[2] * bCurrentRay->Direction[2] - block_min[2];

		TraverseVoxels(block_min,
			       block_x,
			       block_flow_field[i + j + k],
			       t_max[2],
			       bCurrentRay,
			       xyz_Is_1);
	      }

	      t_max[2] = xyz_Is_1[2]
		? t_max[2] + t_delta[2]
		: t_max[2] - t_delta[2];
	    }
	  }
	}
      }

      void RayTracer::BuildClusters()
      {
	ClusterBuilder clusterBuilder(mLatDat, cluster_voxel, cluster_flow_field);
	clusterBuilder.BuildClusters();
	
	mClusters = clusterBuilder.GetClusters();
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
	BuildClusters();
      }

      void RayTracer::Render()
      {
	const float* const projectedUnitX = mScreen->GetUnitVectorProjectionX();
	const float* const projectedUnitY = mScreen->GetUnitVectorProjectionY();

	const float* viewpointCentre = mViewpoint->GetViewpointCentre();

	for (unsigned int clusterId = 0; clusterId < mClusters.size(); clusterId++)
	{
	  const Cluster* thisCluster = &mClusters[clusterId];

	  // the image-based projection of the mClusters bounding box is
	  // calculated here
	  float cluster_x[3];
	  cluster_x[0] = thisCluster->minBlock.x - viewpointCentre[0];
	  cluster_x[1] = thisCluster->minBlock.y - viewpointCentre[1];
	  cluster_x[2] = thisCluster->minBlock.z - viewpointCentre[2];
		    

	  float **block_flow_field = cluster_flow_field[clusterId];

	  float subimageMins[2], subimageMaxes[2];

	  subimageMins[0] = subimageMins[1] = std::numeric_limits<float>::max();
	  subimageMaxes[0] = subimageMaxes[1] = std::numeric_limits<float>::min();

	  float p1[3];

	  // Temp fix due to refactoring of cluster builder

	  float lMinMax_x[2];
	  lMinMax_x[0] = thisCluster->minSite.x;
	  lMinMax_x[1] = thisCluster->maxSite.x;

	  float lMinMax_y[2];
	  lMinMax_y[0] = thisCluster->minSite.y;
	  lMinMax_y[1] = thisCluster->maxSite.y;

	  float lMinMax_z[2];
	  lMinMax_z[0] = thisCluster->minSite.z;
	  lMinMax_z[1] = thisCluster->maxSite.z;  

	  for (int i = 0; i < 2; i++)
	  {
	    p1[0] = lMinMax_x[i];

	    for (int j = 0; j < 2; j++)
	    {
	      p1[1] = lMinMax_y[j];

	      for (int k = 0; k < 2; k++)
	      {
		p1[2] = lMinMax_z[k];

		float p2[3];
		mViewpoint->Project(p1, p2);

		subimageMins[0] = fminf(subimageMins[0], p2[0]);
		subimageMaxes[0] = fmaxf(subimageMaxes[0], p2[0]);

		subimageMins[1] = fminf(subimageMins[1], p2[1]);
		subimageMaxes[1] = fmaxf(subimageMaxes[1], p2[1]);
	      }
	    }
	  }

	  int subimageMinXY[2], subimageMaxXY[2];

	  mScreen->Transform<int> (subimageMins, subimageMinXY);
	  mScreen->Transform<int> (subimageMaxes, subimageMaxXY);

	  // If the entire sub-image is off the screen, continue to the next cluster.
	  if (subimageMinXY[0] >= mScreen->GetPixelsX() || subimageMaxXY[0] < 0 || subimageMinXY[1]
	      >= mScreen->GetPixelsY() || subimageMaxXY[1] < 0)
	  {
	    continue;
	  }

	  // Crop the sub-image to the screen.
	  subimageMinXY[0] = util::NumericalFunctions::max(subimageMinXY[0], 0);
	  subimageMaxXY[0] = util::NumericalFunctions::min(subimageMaxXY[0], mScreen->GetPixelsX()
							   - 1);
	  subimageMinXY[1] = util::NumericalFunctions::max(subimageMinXY[1], 0);
	  subimageMaxXY[1] = util::NumericalFunctions::min(subimageMaxXY[1], mScreen->GetPixelsY()
							   - 1);

	  AABB aabb;
	  aabb.acc_1 = thisCluster->maxSite.x - viewpointCentre[0];
	  aabb.acc_2 = thisCluster->minSite.x - viewpointCentre[0];
	  aabb.acc_3 = thisCluster->maxSite.y - viewpointCentre[1];
	  aabb.acc_4 = thisCluster->minSite.y - viewpointCentre[1];
	  aabb.acc_5 = thisCluster->maxSite.z - viewpointCentre[2];
	  aabb.acc_6 = thisCluster->minSite.z - viewpointCentre[2];

	  float par3[3];
	  const float* vtx = mScreen->GetVtx();
	  for (int l = 0; l < 3; l++)
	  {
	    par3[l] = vtx[l] + (float) subimageMinXY[0] * projectedUnitX[l]
	      + (float) subimageMinXY[1] * projectedUnitY[l];
	  }

	  for (int subImageX = subimageMinXY[0]; subImageX <= subimageMaxXY[0]; ++subImageX)
	  {
	    float lRayDirection[3];
	    for (int l = 0; l < 3; l++)
	    {
	      lRayDirection[l] = par3[l];
	    }

	    for (int subImageY = subimageMinXY[1]; subImageY <= subimageMaxXY[1]; ++subImageY)
	    {
	      Ray lRay;

	      lRay.Direction[0] = lRayDirection[0];
	      lRay.Direction[1] = lRayDirection[1];
	      lRay.Direction[2] = lRayDirection[2];

	      float lInverseDirectionMagnitude = 1.0F / sqrtf(lRayDirection[0] * lRayDirection[0]
							      + lRayDirection[1] * lRayDirection[1] + lRayDirection[2] * lRayDirection[2]);

	      lRay.Direction[0] *= lInverseDirectionMagnitude;
	      lRay.Direction[1] *= lInverseDirectionMagnitude;
	      lRay.Direction[2] *= lInverseDirectionMagnitude;
	      // 
	      lRay.InverseDirection[0] = 1.0F / lRay.Direction[0];
	      lRay.InverseDirection[1] = 1.0F / lRay.Direction[1];
	      lRay.InverseDirection[2] = 1.0F / lRay.Direction[2];

	      bool lRayInPositiveDirection[3];
	      lRayInPositiveDirection[0] = lRay.Direction[0] > 0.0F;
	      lRayInPositiveDirection[1] = lRay.Direction[1] > 0.0F;
	      lRayInPositiveDirection[2] = lRay.Direction[2] > 0.0F;

	      lRayDirection[0] += projectedUnitY[0];
	      lRayDirection[1] += projectedUnitY[1];
	      lRayDirection[2] += projectedUnitY[2];

	      float t_near, t_far;
	      AABBvsRay(&aabb, lRay.InverseDirection, lRayInPositiveDirection, &t_near, &t_far);

	      if (t_near > t_far)
	      {
		continue;
	      }

	      float ray_dx[3];
	      ray_dx[0] = t_near * lRay.Direction[0] - cluster_x[0];
	      ray_dx[1] = t_near * lRay.Direction[1] - cluster_x[1];
	      ray_dx[2] = t_near * lRay.Direction[2] - cluster_x[2];

	      lRay.VelocityColour[0] = 0.0F;
	      lRay.VelocityColour[1] = 0.0F;
	      lRay.VelocityColour[2] = 0.0F;

	      lRay.StressColour[0] = 0.0F;
	      lRay.StressColour[1] = 0.0F;
	      lRay.StressColour[2] = 0.0F;

	      lRay.Length = 0.0F;
	      lRay.MinT = std::numeric_limits<float>::max();
	      lRay.Density = -1.0F;

	      TraverseBlocks(thisCluster, lRayInPositiveDirection, ray_dx, block_flow_field, &lRay);

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
	    par3[0] += projectedUnitX[0];
	    par3[1] += projectedUnitX[1];
	    par3[2] += projectedUnitX[2];
	  }
	}
      }

      void RayTracer::UpdateClusterVoxel(site_t i,
					 distribn_t density,
					 distribn_t velocity,
					 distribn_t stress)
      {
	assert(static_cast<site_t>(static_cast<unsigned int>(i)) == i);

	if(cluster_voxel[3 * i] == NULL)
	{
	  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>
	    ("Failed to access cluster voxel at site %u",(unsigned int) i); 
	  assert(false);
	}

	cluster_voxel[3 * i][0] = (float) density;
	cluster_voxel[3 * i][1] = (float) velocity;
	cluster_voxel[3 * i][2] = (float) stress;
      }

      RayTracer::~RayTracer()
      {
	for (unsigned int n = 0; n < mClusters.size(); n++)
	{
	  for (int m = 0; m < (mClusters[n].blocksX * mClusters[n].blocksY * mClusters[n].blocksZ); m++)
	  {
	    if (cluster_flow_field[n][m] != NULL)
	    {
	      delete[] cluster_flow_field[n][m];
	    }
	  }
	  delete[] cluster_flow_field[n];
	}

	delete[] cluster_flow_field;
	delete[] cluster_voxel;
      }
    }
  }
}
