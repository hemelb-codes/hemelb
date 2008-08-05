#include "config.h"


void (*rtAABBvsRay[2][2][2]) (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far);

void (*rtTraverseBlocks[2][2][2]) (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]));


void rtAABBvsRay000 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}

void rtAABBvsRay001 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}

void rtAABBvsRay010 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}

void rtAABBvsRay011 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}

void rtAABBvsRay100 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}

void rtAABBvsRay101 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}

void rtAABBvsRay110 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}

void rtAABBvsRay111 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
  *t_far  = fminf(tx1, fminf(ty1, tz1));
}


void rtUpdateColour (float dt, float palette[], float col[])
{
  col[0] += dt * palette[0];
  col[1] += dt * palette[1];
  col[2] += dt * palette[2];
}


void rtUpdateRayData (float *flow_field, float ray_t, float ray_segment, void (*ColourPalette) (float value, float col[]))
{
  if (*flow_field < 0.F) return; // solid voxel
  
  float palette[3];
  
  // update the volume rendering of the velocity flow field
  float scaled_velocity = *(flow_field+1) * vis_velocity_threshold_max_inv;
  
  ColourPalette (scaled_velocity, palette);
  
  ray_vel_col[0] += ray_segment * palette[0];
  ray_vel_col[1] += ray_segment * palette[1];
  ray_vel_col[2] += ray_segment * palette[2];
  
  // update the volume rendering of the von Mises stress flow field
  float scaled_stress = *(flow_field+2) * vis_stress_threshold_max_inv;
  
  ColourPalette (scaled_stress, palette);
  
  ray_stress_col[0] += ray_segment * palette[0];
  ray_stress_col[1] += ray_segment * palette[1];
  ray_stress_col[2] += ray_segment * palette[2];
  
  if (ray_density >= 0.F) return;
  
  ray_t_min = ray_t;
  
  // keep track of the density nearest to the view point
  ray_density = *flow_field;
  
  // keep track of the von Mises stress nearest to the view point
  ray_stress = *(flow_field+2);
}


void rtTraverseVoxels000 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}


void rtTraverseVoxels001 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[2] += ray_inv[2];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}


void rtTraverseVoxels010 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[1] += ray_inv[1];
 
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}


void rtTraverseVoxels011 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[0] -= ray_inv[0];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}


void rtTraverseVoxels100 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[0] += ray_inv[0];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}


void rtTraverseVoxels101 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[1] -= ray_inv[1];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}


void rtTraverseVoxels110 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[2] -= ray_inv[2];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}


void rtTraverseVoxels111 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			  void (*ColourPalette) (float value, float col[]))
{
  float t_max[3];
  
  int i_vec[3];
  int i, j, k;
  
  
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < 3; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < 3; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[0]-t, ColourPalette);
	      
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[1]-t, ColourPalette);
	      
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      rtUpdateRayData (&voxel_flow_field[ 3*(i+j+k) ], t, t_max[2]-t, ColourPalette);
	      
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}


void rtTraverseBlocks000 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels000 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels000 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels000 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels000 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels000 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverseBlocks001 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels001 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  t_max[2] += t_delta[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels001 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels001 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels001 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels001 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

void rtTraverseBlocks010 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels010 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  t_max[1] += t_delta[1];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels010 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels010 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels010 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels010 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverseBlocks011 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels011 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  t_max[0] -= t_delta[0];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels011 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels011 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels011 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels011 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

void rtTraverseBlocks100 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels100 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  t_max[0] += t_delta[0];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels100 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels100 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels100 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels100 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverseBlocks101 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels101 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  t_max[1] -= t_delta[1];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels101 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels101 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels101 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels101 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

void rtTraverseBlocks110 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels110 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  t_max[2] -= t_delta[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels110 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels110 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels110 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels110 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverseBlocks111 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[3];
  float t_max[3];
  float block_x[3];
  float t_delta[3];
  float dx[3];
  
  int i_vec[3];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
      
      rtTraverseVoxels111 (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette);
    }
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels111 (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette);
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels111 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels111 (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette);
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
		  
		  rtTraverseVoxels111 (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette);
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}


void rtBuildClusters (Net *net)
{
  int n_x[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, +0, +0, +0, +0, +0, +0, +0, +0, +1, +1, +1, +1, +1, +1, +1, +1, +1};
  int n_y[] = {-1, -1, -1, +0, +0, +0, +1, +1, +1, -1, -1, -1, +0, +0, +1, +1, +1, -1, -1, -1, +0, +0, +0, +1, +1, +1};  
  int n_z[] = {-1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1};
  
  int neigh_i, neigh_j, neigh_k;
  int voxel_min[3], voxel_max[3];
  int block_coord[3], ii[3];
  int cluster_id;
  int block_id;
  int blocks_a, blocks_b;
  int index_a;
  int blocks_buffer_size;
  int are_blocks_incrementing;
  int is_site_found;
  int i, j, k;
  int l, m, n;
  int clusters_max;
  int block_min_x, block_min_y, block_min_z;
  int block_max_x, block_max_y, block_max_z;
  int dummy = 0;
  
  unsigned int my_site_id;
  
  unsigned short int cluster_block_max_i, cluster_block_max_j, cluster_block_max_k;
  
  bool *is_block_visited;
  
  BlockLocation *block_location_a, *block_location_b;
  BlockLocation *block_location_a_p, *block_location_b_p;
  
  DataBlock *map_block_p;
  ProcBlock *proc_block_p;
  
  Cluster *cluster_p;
  
  
  net->cluster_id = (short int *)malloc(sizeof(short int) * blocks);
  
  clusters_max = 20;
  clusters = 0;
  cluster = (Cluster *)malloc(sizeof(Cluster) * clusters_max);
  
  for (n = 0; n < blocks; n++)
    {
      net->cluster_id[ n ] = -1;
    }
  cluster_p = NULL;
  cluster_block_max_i = dummy;
  cluster_block_max_j = dummy;
  cluster_block_max_k = dummy;
  
  is_block_visited = (bool *)malloc(sizeof(bool) * blocks);
  
  blocks_buffer_size = 10000;
  block_location_a = (BlockLocation *)malloc(sizeof(BlockLocation) * blocks_buffer_size);
  block_location_b = (BlockLocation *)malloc(sizeof(BlockLocation) * blocks_buffer_size);
  
  for (n = 0; n < blocks; n++)
    {
      is_block_visited[ n ] = 0;
    }
  block_min_x = +1000000000;
  block_min_y = +1000000000;
  block_min_z = +1000000000;
  block_max_x = -1000000000;
  block_max_y = -1000000000;
  block_max_z = -1000000000;
  
  n = -1;
  
  for (i = 0; i < blocks_x; i++)
    {
      for (j = 0; j < blocks_y; j++)
	{
	  for (k = 0; k < blocks_z; k++)
	    {
	      if ((proc_block_p = &net->proc_block[ ++n ])->proc_id == NULL) continue;
	      
	      block_min_x = min(block_min_x, i);
	      block_min_y = min(block_min_y, j);
	      block_min_z = min(block_min_z, k);
	      block_max_x = max(block_max_x, i);
	      block_max_y = max(block_max_y, j);
	      block_max_z = max(block_max_z, k);
	      
	      if (is_block_visited[ n ]) continue;
	      
	      is_block_visited[ n ] = 1;
	      
	      blocks_a = 0;
	      
	      for (m = 0; m < sites_in_a_block; m++)
		{
		  if (proc_block_p->proc_id[ m ] == net->id)
		    {
		      block_location_a_p = &block_location_a[ 0 ];
		      block_location_a_p->i = i;
		      block_location_a_p->j = j;
		      block_location_a_p->k = k;
		      blocks_a = 1;
		      break;
		    }
		}
	      if (blocks_a == 0) continue;
	      
	      if (clusters == clusters_max)
		{
		  clusters_max <<= 1;
		  cluster = (Cluster *)realloc(cluster, sizeof(Cluster) * clusters_max);
		}
	      cluster_p = &cluster[ clusters ];
	      ++clusters;
	      
	      cluster_p->block_min[0] = i;
	      cluster_p->block_min[1] = j;
	      cluster_p->block_min[2] = k;
	      
	      cluster_block_max_i = i;
	      cluster_block_max_j = j;
	      cluster_block_max_k = k;
	      
	      net->cluster_id[ n ] = clusters - 1;
	      
	      are_blocks_incrementing = 1;
	      
	      while (are_blocks_incrementing)
		{
		  blocks_b = 0;
		  are_blocks_incrementing = 0;
		  
		  for (index_a = 0; index_a < blocks_a; index_a++)
		    {
		      block_location_a_p = &block_location_a[ index_a ];
		      
		      for (l = 0; l < 26; l++)
			{
			  neigh_i = block_location_a_p->i + n_x[ l ];
			  neigh_j = block_location_a_p->j + n_y[ l ];
			  neigh_k = block_location_a_p->k + n_z[ l ];
			  
			  if (neigh_i == -1 || neigh_i == blocks_x) continue;
			  if (neigh_j == -1 || neigh_j == blocks_y) continue;
			  if (neigh_k == -1 || neigh_k == blocks_z) continue;
			  
			  block_id = (neigh_i * blocks_y + neigh_j) * blocks_z + neigh_k;
			  
			  if (is_block_visited[ block_id ] ||
			      (proc_block_p = &net->proc_block[ block_id ])->proc_id == NULL)
			    {
			      continue;
			    }
			  is_site_found = 0;
			  
			  for (m = 0; m < sites_in_a_block; m++)
			    {
			      if (proc_block_p->proc_id[ m ] == net->id)
				{
				  is_site_found = 1;
				  break;
				}
			    }
			  if (!is_site_found) continue;
			  
			  is_block_visited[ block_id ] = 1;
			  
			  are_blocks_incrementing = 1;
			  
			  if (blocks_b == blocks_buffer_size)
			    {
			      blocks_buffer_size *= 2;
			      block_location_a = (BlockLocation *)realloc(block_location_a,
									  sizeof(BlockLocation) * blocks_buffer_size);
			      block_location_b = (BlockLocation *)realloc(block_location_b,
									  sizeof(BlockLocation) * blocks_buffer_size); 
			    }
			  block_location_b_p = &block_location_b[ blocks_b ];
			  block_location_b_p->i = neigh_i;
			  block_location_b_p->j = neigh_j;
			  block_location_b_p->k = neigh_k;
			  ++blocks_b;
			  
			  cluster_p->block_min[0] = min((int)cluster_p->block_min[0], neigh_i);
			  cluster_p->block_min[1] = min((int)cluster_p->block_min[1], neigh_j);
			  cluster_p->block_min[2] = min((int)cluster_p->block_min[2], neigh_k);
			  
			  cluster_block_max_i = max((int)cluster_block_max_i, neigh_i);
			  cluster_block_max_j = max((int)cluster_block_max_j, neigh_j);
			  cluster_block_max_k = max((int)cluster_block_max_k, neigh_k);
			  
			  net->cluster_id[ block_id ] = clusters - 1;
			}
		    }
		  block_location_a_p = block_location_a;
		  block_location_a = block_location_b;
		  block_location_b = block_location_a_p;
		  blocks_a = blocks_b;
		}
	      cluster_p->x[0] = cluster_p->block_min[0] * block_size - 0.5F * sites_x;
	      cluster_p->x[1] = cluster_p->block_min[1] * block_size - 0.5F * sites_y;
	      cluster_p->x[2] = cluster_p->block_min[2] * block_size - 0.5F * sites_z;
	      
	      cluster_p->blocks_x = 1 + cluster_block_max_i - cluster_p->block_min[0];
	      cluster_p->blocks_y = 1 + cluster_block_max_j - cluster_p->block_min[1];
	      cluster_p->blocks_z = 1 + cluster_block_max_k - cluster_p->block_min[2];
	    }
	}
    }
  free(block_location_b);
  free(block_location_a);
  
  free(is_block_visited);
  
  vis_ctr_x = 0.5F * block_size * (block_min_x + block_max_x) - vis->half_dim[0];
  vis_ctr_y = 0.5F * block_size * (block_min_y + block_max_y) - vis->half_dim[1];
  vis_ctr_z = 0.5F * block_size * (block_min_z + block_max_z) - vis->half_dim[2];
  
  
  cluster_voxel = (float **)malloc(sizeof(float *) * 3 * net->my_sites);
  
  cluster_flow_field = (float ***)malloc(sizeof(float **) * clusters);
  
  for (cluster_id = 0; cluster_id < clusters; cluster_id++)
    {
      cluster_p = &cluster[ cluster_id ];
      
      cluster_flow_field[ cluster_id ] = (float **)malloc(sizeof(float *) *
							  cluster_p->blocks_x *
							  cluster_p->blocks_y *
							  cluster_p->blocks_z);
      for (l = 0; l < 3; l++)
	{
	  voxel_min[l] = +1000000000;
	  voxel_max[l] = -1000000000;
	}
      n = -1;
      
      for (i = 0; i < cluster_p->blocks_x; i++)
	{
	  block_coord[0] = (i + cluster_p->block_min[0]) * block_size;
	  
	  for (j = 0; j < cluster_p->blocks_y; j++)
	    {
	      block_coord[1] = (j + cluster_p->block_min[1]) * block_size;
	      
	      for (k = 0; k < cluster_p->blocks_z; k++)
		{
		  block_coord[2] = (k + cluster_p->block_min[2]) * block_size;
		  
		  block_id = ((i + cluster_p->block_min[0]) * blocks_y +
			      (j + cluster_p->block_min[1])) * blocks_z +
		    (k + cluster_p->block_min[2]);
		  
		  cluster_flow_field[ cluster_id ][ ++n ] = NULL;
		  
		  if (net->cluster_id[ block_id ] != cluster_id)
		    {
		      continue;
		    }
		  map_block_p = &net->map_block[ block_id ];
		  
		  cluster_flow_field[ cluster_id ][n] = (float *)malloc(sizeof(float) * 3 * sites_in_a_block);
		  
		  m = -1;
		  
		  for (ii[0] = 0; ii[0] < block_size; ii[0]++)
		    {
		      for (ii[1] = 0; ii[1] < block_size; ii[1]++)
			{
			  for (ii[2] = 0; ii[2] < block_size; ii[2]++)
			    {
			      my_site_id = map_block_p->site_data[ ++m ];
			      
			      if (my_site_id & (1U << 31U))
				{
				  cluster_flow_field[ cluster_id ][n][ 3*m+0 ] = -1.F;
				  cluster_flow_field[ cluster_id ][n][ 3*m+1 ] = -1.F;
				  cluster_flow_field[ cluster_id ][n][ 3*m+2 ] = -1.F;
				  continue;
				}
			      cluster_flow_field[ cluster_id ][n][ 3*m+0 ] = 1.F;
			      cluster_flow_field[ cluster_id ][n][ 3*m+1 ] = 1.F;
			      cluster_flow_field[ cluster_id ][n][ 3*m+2 ] = 1.F;
			      cluster_voxel[ 3*my_site_id+0 ] = &cluster_flow_field[ cluster_id ][n][ 3*m+0 ];
			      cluster_voxel[ 3*my_site_id+1 ] = &cluster_flow_field[ cluster_id ][n][ 3*m+1 ];
			      cluster_voxel[ 3*my_site_id+2 ] = &cluster_flow_field[ cluster_id ][n][ 3*m+2 ];
				  
			      for (l = 0; l < 3; l++)
				{
				  voxel_min[l] = min(voxel_min[l], ii[l] + block_coord[l]);
				  voxel_max[l] = max(voxel_max[l], ii[l] + block_coord[l]);
				}
			    }
			}
		    }
		}
	    }
	}
      cluster_p->minmax_x[0] = (float)voxel_min[0] - 0.5F * (float)sites_x;
      cluster_p->minmax_y[0] = (float)voxel_min[1] - 0.5F * (float)sites_y;
      cluster_p->minmax_z[0] = (float)voxel_min[2] - 0.5F * (float)sites_z;
      
      cluster_p->minmax_x[1] = (float)(voxel_max[0] + 1) - 0.5F * (float)sites_x;
      cluster_p->minmax_y[1] = (float)(voxel_max[1] + 1) - 0.5F * (float)sites_y;
      cluster_p->minmax_z[1] = (float)(voxel_max[2] + 1) - 0.5F * (float)sites_z;
    }
  free(net->cluster_id);
}


void rtInit (Net *net)
{
  rtAABBvsRay[0][0][0] = rtAABBvsRay000;
  rtAABBvsRay[0][0][1] = rtAABBvsRay001;
  rtAABBvsRay[0][1][0] = rtAABBvsRay010;
  rtAABBvsRay[0][1][1] = rtAABBvsRay011;
  rtAABBvsRay[1][0][0] = rtAABBvsRay100;
  rtAABBvsRay[1][0][1] = rtAABBvsRay101;
  rtAABBvsRay[1][1][0] = rtAABBvsRay110;
  rtAABBvsRay[1][1][1] = rtAABBvsRay111;
  
  rtTraverseBlocks[0][0][0] = rtTraverseBlocks000;
  rtTraverseBlocks[0][0][1] = rtTraverseBlocks001;
  rtTraverseBlocks[0][1][0] = rtTraverseBlocks010;
  rtTraverseBlocks[0][1][1] = rtTraverseBlocks011;
  rtTraverseBlocks[1][0][0] = rtTraverseBlocks100;
  rtTraverseBlocks[1][0][1] = rtTraverseBlocks101;
  rtTraverseBlocks[1][1][0] = rtTraverseBlocks110;
  rtTraverseBlocks[1][1][1] = rtTraverseBlocks111;
  
  rtBuildClusters (net);
}


void rtRayTracing (void (*ColourPalette) (float value, float col[]))
{
  // the volume rendering is performed here
  
  float screen_max[4];
  float screen_vtx[4];
  float p0[3], p1[3], p2[3];
  float ray_dx[3];
  float cluster_x[3];
  float dir[3];
  float par1[3], par2[3], par3[3];
  float subimage_vtx[4];
  float scale_vec[4];
  float t_near, t_far;
  float **block_flow_field;
  float temp1;
  
  int subimage_pix[4];
  int ray_sign[3];

  int pixels_x, pixels_y;
  int cluster_id;
  // int viewpoint_flag;
  int i, j, k, l;
  
  AABB aabb;
  
  Cluster *cluster_p;
  
  ColPixel col_pixel;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  screen_max[0] = screen.max_x;
  screen_max[1] = screen.max_x;
  screen_max[2] = screen.max_y;
  screen_max[3] = screen.max_y;
  
  for (l = 0; l < 3; l++)
    {
      p0[l] = viewpoint.x[l];
    }
  for (l = 0; l < 3; l++)
    {
      par1[l] = screen.dir1[l];
      par2[l] = screen.dir2[l];
      screen_vtx[l] = screen.vtx[l] + 0.5F * par1[l] + 0.5F * par2[l];
    }
  scale_vec[0] = scale_vec[1] = screen.scale_x;
  scale_vec[2] = scale_vec[3] = screen.scale_y;
  
  for (cluster_id = 0; cluster_id < clusters; cluster_id++)
    {
      cluster_p = &cluster[ cluster_id ];
      
      // the image-based projection of the cluster bounding box is
      // calculated here
      
      for (l = 0; l < 3; l++)
	{
	  cluster_x[l] = cluster_p->x[l] - p0[l];
	}
      cluster_blocks_vec[0] = cluster_p->blocks_x - 1;
      cluster_blocks_vec[1] = cluster_p->blocks_y - 1;
      cluster_blocks_vec[2] = cluster_p->blocks_z - 1;
      cluster_blocks_z      = cluster_p->blocks_z;
      cluster_blocks_yz     = (int)cluster_p->blocks_y * (int)cluster_p->blocks_z;
      cluster_blocks        = (int)cluster_p->blocks_x * cluster_blocks_yz;
      
      block_flow_field = cluster_flow_field[ cluster_id ];
      
      subimage_vtx[0] =  1.e+30F;
      subimage_vtx[1] = -1.e+30F;
      subimage_vtx[2] =  1.e+30F;
      subimage_vtx[3] = -1.e+30F;
      
      for (i = 0; i < 2; i++)
	{
	  p1[0] = cluster_p->minmax_x[i];
	  
	  for (j = 0; j < 2; j++)
	    {
	      p1[1] = cluster_p->minmax_y[j];
	      
	      for (k = 0; k < 2; k++)
		{
		  p1[2] = cluster_p->minmax_z[k];
		  
		  visProject (p1, p2);
		  
		  subimage_vtx[0] = fminf(subimage_vtx[0], p2[0]);
		  subimage_vtx[1] = fmaxf(subimage_vtx[1], p2[0]);
		  subimage_vtx[2] = fminf(subimage_vtx[2], p2[1]);
		  subimage_vtx[3] = fmaxf(subimage_vtx[3], p2[1]);
		}
	    }
	}
      subimage_pix[0] = (int)(scale_vec[0] * (subimage_vtx[0] + screen_max[0]));
      subimage_pix[1] = (int)(scale_vec[1] * (subimage_vtx[1] + screen_max[1]));
      subimage_pix[2] = (int)(scale_vec[2] * (subimage_vtx[2] + screen_max[2]));
      subimage_pix[3] = (int)(scale_vec[3] * (subimage_vtx[3] + screen_max[3]));
      
      if (subimage_pix[0] >= pixels_x || subimage_pix[1] < 0 ||
	  subimage_pix[2] >= pixels_y || subimage_pix[3] < 0)
	{
	  continue;
	}
      subimage_pix[0] = max(subimage_pix[0], 0);
      subimage_pix[1] = min(subimage_pix[1], pixels_x - 1);
      subimage_pix[2] = max(subimage_pix[2], 0);
      subimage_pix[3] = min(subimage_pix[3], pixels_y - 1);
      
      // if (p0[0] >= v[0][0] && p0[1] >= v[0][1] && p0[2] >= v[0][2] &&
      // 	  p0[0] <= v[1][0] && p0[1] <= v[1][1] && p0[2] <= v[1][2])
      // 	{
      // 	  viewpoint_flag = 1;
      // 	}
      // else
      // 	{
      // 	  viewpoint_flag = 0;
      // 	}
      aabb.acc_1 = cluster_p->minmax_x[1] - p0[0];
      aabb.acc_2 = cluster_p->minmax_x[0] - p0[0];
      aabb.acc_3 = cluster_p->minmax_y[1] - p0[1];
      aabb.acc_4 = cluster_p->minmax_y[0] - p0[1];
      aabb.acc_5 = cluster_p->minmax_z[1] - p0[2];
      aabb.acc_6 = cluster_p->minmax_z[0] - p0[2];
      
      for (l = 0; l < 3; l++)
	{
	  par3[l] = screen_vtx[l] + subimage_pix[0] * par1[l] + subimage_pix[2] * par2[l];
	}
      for (i = subimage_pix[0]; i <= subimage_pix[1]; i++)
	{
	  for (l = 0; l < 3; l++)
	    {
	      dir[l] = par3[l];
	    }
	  for (j = subimage_pix[2]; j <= subimage_pix[3]; j++)
	    {
	      ray_dir[0] = dir[0];
	      ray_dir[1] = dir[1];
	      ray_dir[2] = dir[2];
	      
	      temp1 = 1.F / sqrtf(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
	      
	      ray_dir[0] *= temp1;
	      ray_dir[1] *= temp1;
	      ray_dir[2] *= temp1;
	      
	      ray_inv[0] = 1.F / ray_dir[0];
	      ray_inv[1] = 1.F / ray_dir[1];
	      ray_inv[2] = 1.F / ray_dir[2];
	      
	      ray_sign[0] = ray_dir[0] > 0.F;
	      ray_sign[1] = ray_dir[1] > 0.F;
	      ray_sign[2] = ray_dir[2] > 0.F;
	      
	      dir[0] += par2[0];
	      dir[1] += par2[1];
	      dir[2] += par2[2];
	      
	      // t_near = 0.F;
	      
	      // if (!viewpoint_flag)
		{
		  (*rtAABBvsRay[ray_sign[0]][ray_sign[1]][ray_sign[2]])
		    (&aabb, ray_inv[0], ray_inv[1], ray_inv[2], &t_near, &t_far);
		  
		  if (t_near > t_far) continue;
		  
		  ray_dx[0] = t_near * ray_dir[0] - cluster_x[0];
		  ray_dx[1] = t_near * ray_dir[1] - cluster_x[1];
		  ray_dx[2] = t_near * ray_dir[2] - cluster_x[2];
		}
	      // else
	      // 	{
	      // 	  ray_dx[0] = cluster_x[0];
	      // 	  ray_dx[1] = cluster_x[1];
	      // 	  ray_dx[2] = cluster_x[2];
	      // 	}
	      ray_vel_col[0] = 0.F;
	      ray_vel_col[1] = 0.F;
	      ray_vel_col[2] = 0.F;
	      ray_stress_col[0] = 0.F;
	      ray_stress_col[1] = 0.F;
	      ray_stress_col[2] = 0.F;
	      ray_t_min = 1.e+30F;
	      ray_density = -1.F;
	      
	      (*rtTraverseBlocks[ray_sign[0]][ray_sign[1]][ray_sign[2]])
		(ray_dx, block_flow_field, ColourPalette);
	      
	      if (ray_t_min >= 1.e+30F) continue;
	      
	      col_pixel.vel_r    = ray_vel_col[ 0 ];
	      col_pixel.vel_g    = ray_vel_col[ 1 ];
	      col_pixel.vel_b    = ray_vel_col[ 2 ];
	      col_pixel.stress_r = ray_stress_col[ 0 ];
	      col_pixel.stress_g = ray_stress_col[ 1 ];
	      col_pixel.stress_b = ray_stress_col[ 2 ];
	      col_pixel.t        = ray_t_min + t_near;
	      col_pixel.density  = ray_density;
	      col_pixel.stress   = ray_stress;
	      col_pixel.i        = i * (1 << 16) + j;
	      
	      visWritePixel (&col_pixel);
	    }
	  par3[0] += par1[0];
	  par3[1] += par1[1];
	  par3[2] += par1[2];
	}
    }
}


void rtEnd (void)
{
  int m, n;
  
  
  for (n = 0; n < clusters; n++)
    {
      for (m = 0; m < cluster[n].blocks_x * cluster[n].blocks_y * cluster[n].blocks_z; m++)
  	{
  	  if (cluster_flow_field[ n ][ m ] != NULL)
  	    {
  	      free(cluster_flow_field[ n ][ m ]);
  	    }
  	}
      free(cluster_flow_field[ n ]);
    }
  free(cluster_flow_field);
  
  free(cluster_voxel);
  
  free(cluster);
}


void visRotate (float sin_1, float cos_1,
		float sin_2, float cos_2,
		float  x1, float  y1, float  z1,
		float *x2, float *y2, float *z2)
{
  //sin_1 = sin(longitude), cos_1 = cos(longitude)
  //sin_2 = sin(latitude) , cos_2 = cos(latitude)
  //rotation of latitude  DEG around the X axis followed by a
  //rotation of longitude DEG around the Y axis
  
  float temp = z1 * cos_2 - y1 * sin_2;
  
  *x2 = temp * sin_1 + x1 * cos_1;
  *y2 =   z1 * sin_2 + y1 * cos_2;
  *z2 = temp * cos_1 - x1 * sin_1;
}


void visProject (float p1[], float p2[])
{
  float x1[3], x2[3];
  float temp;
  
  int l;
  
  
  for (l = 0; l < 3; l++)
    {
      x1[l] = p1[l] - viewpoint.x[l];
    }
  temp = viewpoint.cos_1 * x1[2] + viewpoint.sin_1 * x1[0];
  
  x2[0] = viewpoint.cos_1 * x1[0] - viewpoint.sin_1 * x1[2];
  x2[1] = viewpoint.cos_2 * x1[1] - viewpoint.sin_2 * temp;
  x2[2] = viewpoint.cos_2 * temp + viewpoint.sin_2 * x1[1];
  
  temp = viewpoint.dist / (-x2[2]);
  
  for (l = 0; l < 3; l++)
    {
      p2[l] = temp * x2[l];
    }
  p2[2] = -x2[2];
}


void visProjection (float ortho_x, float ortho_y,
		    int pixels_x, int pixels_y,
		    float ctr_x, float ctr_y, float ctr_z,
		    float rad,
		    float longitude, float latitude,
		    float dist,
		    float zoom)
{
  float temp;
  
  
  screen.max_x = ortho_x / zoom;
  screen.max_y = ortho_y / zoom;
  
  screen.pixels_x = pixels_x;
  screen.pixels_y = pixels_y;
  
  temp = longitude * 0.01745329F;
  
  viewpoint.sin_1 = sinf(temp);
  viewpoint.cos_1 = cosf(temp);
  
  temp = latitude * 0.01745329F;
  
  viewpoint.sin_2 = sinf(temp);
  viewpoint.cos_2 = cosf(temp);
  
  temp = rad * viewpoint.cos_2;
  
  viewpoint.x[0] = temp * viewpoint.sin_1 + ctr_x;
  viewpoint.x[1] = rad  * viewpoint.sin_2 + ctr_y;
  viewpoint.x[2] = temp * viewpoint.cos_1 + ctr_z;
  
  viewpoint.dist = dist;
  
  temp = dist / rad;
  
  ctr_x = viewpoint.x[0] + temp * (ctr_x - viewpoint.x[0]);
  ctr_y = viewpoint.x[1] + temp * (ctr_y - viewpoint.x[1]);
  ctr_z = viewpoint.x[2] + temp * (ctr_z - viewpoint.x[2]);
  
  screen.zoom = zoom;
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     screen.max_x, 0.0F, 0.0F,
	     &screen.dir1[0], &screen.dir1[1], &screen.dir1[2]);
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     0.0F, screen.max_y, 0.0F,
	     &screen.dir2[0], &screen.dir2[1], &screen.dir2[2]);
  
  screen.scale_x = (float)pixels_x / (2.F * screen.max_x);
  screen.scale_y = (float)pixels_y / (2.F * screen.max_y);
  
  screen.vtx[0] = ctr_x - screen.dir1[0] - screen.dir2[0] - viewpoint.x[0];
  screen.vtx[1] = ctr_y - screen.dir1[1] - screen.dir2[1] - viewpoint.x[1];
  screen.vtx[2] = ctr_z - screen.dir1[2] - screen.dir2[2] - viewpoint.x[2];
  
  screen.dir1[0] *= (2.F / (float)pixels_x);
  screen.dir1[1] *= (2.F / (float)pixels_x);
  screen.dir1[2] *= (2.F / (float)pixels_x);
  
  screen.dir2[0] *= (2.F / (float)pixels_y);
  screen.dir2[1] *= (2.F / (float)pixels_y);
  screen.dir2[2] *= (2.F / (float)pixels_y);
}


void visWritePixel (ColPixel *col_pixel_p)
{
  int *col_pixel_id_p, i, j;
  
  
  i = PixelI(col_pixel_p->i);
  j = PixelJ(col_pixel_p->i);
  
  if (*(col_pixel_id_p = &col_pixel_id[ i * screen.pixels_y + j ]) != -1)
    {
      col_pixel_send[ *col_pixel_id_p ].vel_r += col_pixel_p->vel_r;
      col_pixel_send[ *col_pixel_id_p ].vel_g += col_pixel_p->vel_g;
      col_pixel_send[ *col_pixel_id_p ].vel_b += col_pixel_p->vel_b;
      
      col_pixel_send[ *col_pixel_id_p ].stress_r += col_pixel_p->stress_r;
      col_pixel_send[ *col_pixel_id_p ].stress_g += col_pixel_p->stress_g;
      col_pixel_send[ *col_pixel_id_p ].stress_b += col_pixel_p->stress_b;
      
      if (col_pixel_p->t < col_pixel_send[ *col_pixel_id_p ].t)
	{
	  col_pixel_send[ *col_pixel_id_p ].t       = col_pixel_p->t;
	  col_pixel_send[ *col_pixel_id_p ].density = col_pixel_p->density;
	  col_pixel_send[ *col_pixel_id_p ].stress  = col_pixel_p->stress;
	}
    }
  else
    {
      if (col_pixels >= COLOURED_PIXELS_PER_PROC_MAX)
	{
	  printf (" too many coloured pixels per proc\n");
	  printf (" the execution is terminated\n");
#ifndef NOMPI
	  MPI_Abort (MPI_COMM_WORLD, 1);
#else
	  exit(1);
#endif
	}
      if (col_pixels == col_pixels_max)
	{
	  col_pixels_max *= 2;
	  // col_pixel_send = (ColPixel *)realloc(col_pixel_send,
	  // 				       sizeof(ColPixel) * col_pixels_max * max(1, (net_machines - 1)));
	  col_pixel_recv = (ColPixel *)realloc(col_pixel_recv,
					       sizeof(ColPixel) * col_pixels_max);
	}
      *col_pixel_id_p = col_pixels;
      
      col_pixel_send[ col_pixels ].vel_r = col_pixel_p->vel_r;
      col_pixel_send[ col_pixels ].vel_g = col_pixel_p->vel_g;
      col_pixel_send[ col_pixels ].vel_b = col_pixel_p->vel_b;
      
      col_pixel_send[ col_pixels ].stress_r = col_pixel_p->stress_r;
      col_pixel_send[ col_pixels ].stress_g = col_pixel_p->stress_g;
      col_pixel_send[ col_pixels ].stress_b = col_pixel_p->stress_b;
      
      col_pixel_send[ col_pixels ].t       = col_pixel_p->t;
      col_pixel_send[ col_pixels ].density = col_pixel_p->density;
      col_pixel_send[ col_pixels ].stress  = col_pixel_p->stress;
      
      col_pixel_send[ col_pixels ].i = col_pixel_p->i;
      ++col_pixels;
    }
}


void xdrWritePixel (ColPixel *col_pixel_p, XDR *xdr_p, void (*ColourPalette) (float value, float col[]))
{
  float col[3];
  
  int bits_per_char = sizeof(char) * 8;
  int pixel_i, pixel_j;
  
  unsigned int pixel_id;
  unsigned int col_data[3];
  
  unsigned char r1, g1, b1;
  unsigned char r2, g2, b2;
  unsigned char r3, g3, b3;
  unsigned char r4, g4, b4;
  
  
  // store pixel id
  pixel_i = PixelI (col_pixel_p->i);
  pixel_j = PixelJ (col_pixel_p->i);
  
  pixel_id = (pixel_i << (2*bits_per_char)) + pixel_j;
  
  // store velocity flow field
  r1 = (unsigned char)max(0, min(255, (int)col_pixel_p->vel_r));
  g1 = (unsigned char)max(0, min(255, (int)col_pixel_p->vel_g));
  b1 = (unsigned char)max(0, min(255, (int)col_pixel_p->vel_b));
  
  // store von Mises stress flow field
  r2 = (unsigned char)max(0, min(255, (int)col_pixel_p->stress_r));
  g2 = (unsigned char)max(0, min(255, (int)col_pixel_p->stress_g));
  b2 = (unsigned char)max(0, min(255, (int)col_pixel_p->stress_b));
  
  // store external density flow field at the wall
  ColourPalette (col_pixel_p->density, col);
  
  r3 = (unsigned char)max(0, min(255, (int)(255.F * col[0])));
  g3 = (unsigned char)max(0, min(255, (int)(255.F * col[1])));
  b3 = (unsigned char)max(0, min(255, (int)(255.F * col[2])));
  
  // store von Mises stress flow field at the wall
  ColourPalette (col_pixel_p->stress, col);
  
  r4 = (unsigned char)max(0, min(255, (int)(255.F * col[0])));
  g4 = (unsigned char)max(0, min(255, (int)(255.F * col[1])));
  b4 = (unsigned char)max(0, min(255, (int)(255.F * col[2])));
  
  xdr_u_int (xdr_p, &pixel_id);
  
  col_data[0] = (r1<<(3*bits_per_char)) + (g1<<(2*bits_per_char)) + (b1<<bits_per_char) + r2;
  col_data[1] = (g2<<(3*bits_per_char)) + (b2<<(2*bits_per_char)) + (r3<<bits_per_char) + g3;
  col_data[2] = (b3<<(3*bits_per_char)) + (r4<<(2*bits_per_char)) + (g4<<bits_per_char) + b4;
  
  xdr_u_int (xdr_p, &col_data[0]);
  xdr_u_int (xdr_p, &col_data[1]);
  xdr_u_int (xdr_p, &col_data[2]);
}


void visConvertThresholds (float physical_velocity_max, float physical_stress_max,
			   float *lattice_velocity_max, float *lattice_stress_max,
			   LBM *lbm)
{
  float physical_kinematic_viscosity = BLOOD_VISCOSITY / BLOOD_DENSITY;
  
  float useful_factor = ((lbm->tau - 0.5) / 3.) / physical_kinematic_viscosity;
  
  
  *lattice_velocity_max = physical_velocity_max * useful_factor * lbm->voxel_size;
  
  *lattice_stress_max = physical_stress_max * (useful_factor * useful_factor *
					       lbm->voxel_size * lbm->voxel_size) / BLOOD_DENSITY;
}


void visReadParameters (char *parameters_file_name, LBM *lbm, Net *net, Vis *vis)
{
  FILE *parameters_file;
  
  float par_to_send[9];
  float ctr_x, ctr_y, ctr_z;
  float longitude, latitude;
  float zoom;
  float density_min, density_max, velocity_max, stress_max;
  float physical_velocity_max, physical_stress_max;
  
  int i;
  
  
  if (net->id == 0)
    {
      fprintf(stderr, "opening ray tracing configuration file %s\n", parameters_file_name);

      parameters_file = fopen (parameters_file_name, "r");

      if( parameters_file == NULL ) {
        fprintf(stderr, "unable to open file %s, exiting\n", parameters_file_name);
        fflush(0x0);
        exit(0x0);
      } else {
        fprintf(stderr, "done\n");
      }

      fflush(NULL);
      
      fscanf (parameters_file, "%e \n", &ctr_x);
      fscanf (parameters_file, "%e \n", &ctr_y);
      fscanf (parameters_file, "%e \n", &ctr_z);
      fscanf (parameters_file, "%e \n", &longitude);
      fscanf (parameters_file, "%e \n", &latitude);
      fscanf (parameters_file, "%e \n", &zoom);
      fscanf (parameters_file, "%e \n", &vis_brightness);
      fscanf (parameters_file, "%e \n", &physical_velocity_max);
      fscanf (parameters_file, "%e \n", &physical_stress_max);
      
      fclose (parameters_file);
      
      visConvertThresholds (physical_velocity_max, physical_stress_max,
			    &velocity_max, &stress_max, lbm);
      
      par_to_send[ 0 ] = ctr_x;
      par_to_send[ 1 ] = ctr_y;
      par_to_send[ 2 ] = ctr_z;
      par_to_send[ 3 ] = longitude;
      par_to_send[ 4 ] = latitude;
      par_to_send[ 5 ] = zoom;
      par_to_send[ 6 ] = vis_brightness;
      par_to_send[ 7 ] = velocity_max;
      par_to_send[ 8 ] = stress_max;
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 9, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  ctr_x          = par_to_send[ 0 ];
  ctr_y          = par_to_send[ 1 ];
  ctr_z          = par_to_send[ 2 ];
  longitude      = par_to_send[ 3 ];
  latitude       = par_to_send[ 4 ];
  zoom           = par_to_send[ 5 ];
  vis_brightness = par_to_send[ 6 ];
  velocity_max   = par_to_send[ 7 ];
  stress_max     = par_to_send[ 8 ];
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
		 PIXELS_X, PIXELS_Y,
		 ctr_x, ctr_y, ctr_z,
		 5.F * vis->system_size,
		 longitude, latitude,
		 0.5F * (5.F * vis->system_size),
		 zoom);
  
  density_min = +1.e+30F;
  density_max = -1.e+30F;
  
  for (i = 0; i < lbm->inlets; i++)
    {
      density_min = fminf(density_min, inlet_density_avg[ i ] - inlet_density_amp[ i ]);
      density_max = fmaxf(density_max, inlet_density_avg[ i ] + inlet_density_amp[ i ]);
      
    }
  for (i = 0; i < lbm->outlets; i++)
    {
      density_min = fminf(density_min, outlet_density_avg[ i ] - outlet_density_amp[ i ]);
      density_max = fmaxf(density_max, outlet_density_avg[ i ] + outlet_density_amp[ i ]);
    }
  vis_density_threshold_min = density_min;
  vis_density_threshold_minmax_inv = 1.F / (density_max - density_min);
  
  vis_velocity_threshold_max_inv = 1.F / velocity_max;
  vis_stress_threshold_max_inv   = 1.F / stress_max;
  
  vis_image_freq = lbm->period / 100;
}
 

void visInit (Net *net, Vis *vis)
{
  blocks_yz = blocks_y * blocks_z;
  
  vis->half_dim[0] = 0.5F * (float)sites_x;
  vis->half_dim[1] = 0.5F * (float)sites_y;
  vis->half_dim[2] = 0.5F * (float)sites_z;
  
  vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1], vis->half_dim[2]));
  
  block_size_f = (float)block_size;
  block_size2 = block_size * block_size;
  block_size3 = block_size * block_size2;
  block_size_1 = block_size - 1;
  
  block_size_inv = 1.F / (float)block_size;
  
#ifndef NOMPI
  int col_pixel_count = 11;
  int col_pixel_blocklengths[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  MPI_Datatype col_pixel_types[11] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				      MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				      MPI_FLOAT,
				      MPI_FLOAT,
				      MPI_FLOAT,
				      MPI_INT,
				      MPI_UB};
  
  MPI_Aint col_pixel_disps[11];
#endif
  
  col_pixels_max = IMAGE_SIZE / 4;
  
  // col_pixel_send = (ColPixel *)malloc(sizeof(ColPixel) *  col_pixels_max * max(1, (net_machines - 1)));
  col_pixel_recv = (ColPixel *)malloc(sizeof(ColPixel) * col_pixels_max);
  //col_pixel_lock = (ColPixel *)malloc(sizeof(ColPixel) * col_pixels_max);
  
  vis_pixels_max = IMAGE_SIZE;
  col_pixel_id = (int *)malloc(sizeof(int) * vis_pixels_max);
  
  for (int i = 0; i < IMAGE_SIZE; i++)
    {
      col_pixel_id[ i ] = -1;
    }
  
#ifndef NOMPI
  // create the derived datatype for the MPI communications
  
  col_pixel_disps[0] = 0;
  
  for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(float) * col_pixel_blocklengths[i - 1]);
	}
      else if (col_pixel_types[i - 1] == MPI_INT)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
	}
    }
  MPI_Type_struct (col_pixel_count, col_pixel_blocklengths, col_pixel_disps, col_pixel_types, &MPI_col_pixel_type);
  MPI_Type_commit (&MPI_col_pixel_type);
#endif
  
  rtInit (net);
}


void visCompositeImage (Net *net)
{
  // here, intra-machine communications are handled through a binary
  // tree pattern and parallel pairwise blocking communications. The
  // master processor of the current machine gets the sub-images of
  // all the processors of that machine. Inter-machine communications,
  // needed if the number of machines is greater than one, take place
  // in the routine visRenderB.
  
  int *col_pixel_id_p;
  int col_pixels_temp;
  int comm_inc, send_id, recv_id;
  int machine_id, master_proc_id;
  int pixels_y;
  int i, j;
  int m, n;
  
  ColPixel *col_pixel1, *col_pixel2;
  
  
  pixels_y = screen.pixels_y;
  
  memcpy (col_pixel_recv, col_pixel_send, col_pixels * sizeof(ColPixel));
  
  // "master_proc_id" will be the identifier of the processor with
  // lowest rank in its machine
  
  master_proc_id = 0;
  
  for (m = 0; m < net->machine_id[ net->id ]; m++)
    {
      master_proc_id += net->procs_per_machine[ m ];
    }
  comm_inc = 1;
  m = 1;
  
  machine_id = net->machine_id[ net->id ];
  
  while (m < net->procs_per_machine[ machine_id ])
    {
      m <<= 1;
      
      for (recv_id = master_proc_id;
	   recv_id < master_proc_id + net->procs_per_machine[ machine_id ];)
	{
	  send_id = recv_id + comm_inc;
	  
	  if (net->id != recv_id && net->id != send_id)
	    {
	      recv_id += comm_inc << 1;
	      continue;
	    }
	  if (send_id >= master_proc_id + net->procs_per_machine[ machine_id ] ||
	      recv_id == send_id)
	    {
	      recv_id += comm_inc << 1;
	      continue;
	    }
	  if (net->id == send_id)
	    {
#ifndef NOMPI
	      net->err = MPI_Send (&col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
	      if (col_pixels > 0)
		{
#ifndef NOMPI
		  net->err = MPI_Send (&col_pixel_send, col_pixels, MPI_col_pixel_type,
				       recv_id, 20, MPI_COMM_WORLD);
#endif
		}
	    }
	  else
	    {
#ifndef NOMPI
	      net->err = MPI_Recv (&col_pixels_temp, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
				   net->status);
	      if (col_pixels_temp > 0)
		{
		  net->err = MPI_Recv (&col_pixel_send, col_pixels_temp, MPI_col_pixel_type,
				       send_id, 20, MPI_COMM_WORLD, net->status);
		}
#else
	      col_pixels_temp = 0;
#endif
	      for (n = 0; n < col_pixels_temp; n++)
		{
		  col_pixel1 = &col_pixel_send[ n ];
		  i = PixelI (col_pixel1->i);
		  j = PixelJ (col_pixel1->i);
		  
		  if (*(col_pixel_id_p = &col_pixel_id[ i * pixels_y + j ]) == -1)
		    {
		      if (col_pixels == col_pixels_max)
			{
			  col_pixels_max <<= 1;
			  col_pixel_recv = (ColPixel *)realloc(col_pixel_recv,
							       sizeof(ColPixel) * col_pixels_max);
			  //col_pixel_lock = (ColPixel *)realloc(col_pixel_lock,
			  //				       sizeof(ColPixel) * col_pixels_max);
			}
		      col_pixel2 = &col_pixel_recv[ *col_pixel_id_p = col_pixels ];
		      
		      memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		      ++col_pixels;
		    }
		  else
		    {
		      col_pixel2 = &col_pixel_recv[ *col_pixel_id_p ];
		      
		      col_pixel2->vel_r += col_pixel1->vel_r;
		      col_pixel2->vel_g += col_pixel1->vel_g;
		      col_pixel2->vel_b += col_pixel1->vel_b;
		      
		      col_pixel2->stress_r += col_pixel1->stress_r;
		      col_pixel2->stress_g += col_pixel1->stress_g;
		      col_pixel2->stress_b += col_pixel1->stress_b;
		      
		      if (col_pixel1->t < col_pixel2->t)
			{
			  col_pixel2->t       = col_pixel1->t;
			  col_pixel2->density = col_pixel1->density;
			  col_pixel2->stress  = col_pixel1->stress;
			}
		    }
		}
	    }
	  if (m < net->procs_per_machine[ machine_id ])
	    {
	      if (net->id == recv_id)
		{
		  memcpy (col_pixel_send, col_pixel_recv,
			  col_pixels * sizeof(ColPixel));
		}
	    }
	  recv_id += comm_inc << 1;
	}
      comm_inc <<= 1;
    }
}


void visRenderA (void (*ColourPalette) (float value, float col[]), Net *net)
{
  int pixels_x, pixels_y;
  int send_id, recv_id;
  int master_proc_id;
  int m;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  if (pixels_x * pixels_y > vis_pixels_max)
    {
      vis_pixels_max = pixels_x * pixels_y;
      
      col_pixel_id = (int *)realloc(col_pixel_id, sizeof(int) * vis_pixels_max);
    }
  col_pixels = 0;
  
  rtRayTracing (ColourPalette);
  
  if (!vis_compositing)
    {
      for (m = 0; m < col_pixels; m++)
	{
	  col_pixel_id[ (PixelI (col_pixel_send[m].i)) * pixels_y + (PixelJ (col_pixel_send[m].i)) ] = -1;
	}
      return;
    }
  
  visCompositeImage (net);
  
  master_proc_id = 0;
  
  for (m = 0; m < net->machine_id[ net->id ]; m++)
    {
      master_proc_id += net->procs_per_machine[ m ];
    }
  
  if (net_machines == 1 || (net->id != 0 && net->id != master_proc_id))
    {
      return;
    }
  
  // inter-machine communications of sub-images begin here
  
  if (net->id != 0)
    {
      recv_id = 0;
#ifndef NOMPI
      net->err = MPI_Send (&col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
      if (col_pixels > 0)
	{
	  memcpy (col_pixel_send, col_pixel_recv,
		  col_pixels * sizeof(ColPixel));
#ifndef NOMPI
	  net->err = MPI_Isend (col_pixel_send,
				col_pixels, MPI_col_pixel_type,
				recv_id, 30, MPI_COMM_WORLD,
				&net->req[ 1 ][ recv_id ]);
#endif
	}
    }
  else
    {
      send_id = net->procs_per_machine[ net->id ];
      
      for (m = 1; m < net_machines; m++)
	{
#ifndef NOMPI
	  net->err = MPI_Recv (&col_pixels_recv[ m-1 ], 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
			       net->status);
#endif
	  if (col_pixels_recv[ m-1 ] > 0)
	    {
#ifndef NOMPI
	      net->err = MPI_Irecv (&col_pixel_send[ (m-1) * (COLOURED_PIXELS_PER_PROC_MAX * sizeof(ColPixel)) ],
				    col_pixels_recv[ m-1 ], MPI_col_pixel_type,
				    send_id, 30, MPI_COMM_WORLD,
				    &net->req[ 1 ][ net->procs + send_id ]);
#endif
	    }
	  send_id += net->procs_per_machine[ m ];
	}
    }
}


void visRenderB (int write_image, char *image_file_name,
		 void (*ColourPalette) (float value, float col[]), Net *net)
{
  // here, the intra-machine communications take place and the buffer
  // to stream to the client or the output image are set
  
  int pixels_y;
  int i, j;
  int m, n;
  int *col_pixel_id_p;
  int send_id, recv_id;
  int master_proc_id;
  int offset;

  ColPixel *col_pixel1, *col_pixel2;
  
  
  if (!vis_compositing) return;
  
  pixels_y = screen.pixels_y;
  
  if (net_machines > 1)
    {
      master_proc_id = 0;
      
      for (m = 0; m < net->machine_id[ net->id ]; m++)
	{
	  master_proc_id += net->procs_per_machine[ m ];
	}
      if (net->id != 0 && net->id != master_proc_id)
	{
	  for (m = 0; m < IMAGE_SIZE; m++)
	    {
	      col_pixel_id[ m ] = -1;
	    }
	  return;
	}
      
      if (net->id != 0)
	{
	  recv_id = 0;
#ifndef NOMPI
	  net->err = MPI_Wait (&net->req[ 1 ][ recv_id ], net->status);
#endif
	}
      else
	{
	  send_id = net->procs_per_machine[ net->id ];
	  
	  for (m = 1; m < net_machines; m++)
	    {
	      if (col_pixels_recv[ m-1 ] == 0)
		{
		  send_id += net->procs_per_machine[ m ];
		  continue;
		}
#ifndef NOMPI
	      net->err = MPI_Wait (&net->req[ 1 ][ net->procs + send_id ], net->status);
#endif
	      offset = (m-1) * COLOURED_PIXELS_PER_PROC_MAX;
	      
	      for (n = 0; n < col_pixels_recv[ m-1 ]; n++)
		{
		  col_pixel1 = &col_pixel_send[ offset + n ];
		  i = PixelI (col_pixel1->i);
		  j = PixelJ (col_pixel1->i);
		  
		  if (*(col_pixel_id_p = &col_pixel_id[ i * pixels_y + j ]) == -1)
		    {
		      if (col_pixels == col_pixels_max)
			{
			  col_pixels_max <<= 1;
			  col_pixel_recv = (ColPixel *)realloc(col_pixel_recv,
							       sizeof(ColPixel) * col_pixels_max);
			}
		      col_pixel2 = &col_pixel_recv[ *col_pixel_id_p = col_pixels ];
		      
		      memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		      ++col_pixels;
		    }
		  else
		    {
		      col_pixel2 = &col_pixel_recv[ *col_pixel_id_p ];
		      
		      col_pixel2->vel_r += col_pixel1->vel_r;
		      col_pixel2->vel_g += col_pixel1->vel_g;
		      col_pixel2->vel_b += col_pixel1->vel_b;
		      
		      col_pixel2->stress_r += col_pixel1->stress_r;
		      col_pixel2->stress_g += col_pixel1->stress_g;
		      col_pixel2->stress_b += col_pixel1->stress_b;
		      
		      if (col_pixel1->t < col_pixel2->t)
			{
			  col_pixel2->t       = col_pixel1->t;
			  col_pixel2->density = col_pixel1->density;
			  col_pixel2->stress  = col_pixel1->stress;
			}
		    }
		}
	    }
	}
    }
  for (m = 0; m < col_pixels; m++)
    {
      col_pixel_id[ (PixelI (col_pixel_recv[m].i)) * pixels_y + (PixelJ (col_pixel_recv[m].i)) ] = -1;
    }
  
  if (is_bench) return;
  
  if (net->id != 0) return;
  
  FILE *image_file;
  XDR	xdr_image_file;
  
  
  for (n = 0; n < col_pixels; n++)
    {
      col_pixel_recv[ n ].vel_r *= (255.F * vis_brightness);
      col_pixel_recv[ n ].vel_g *= (255.F * vis_brightness);
      col_pixel_recv[ n ].vel_b *= (255.F * vis_brightness);
      
      col_pixel_recv[ n ].stress_r *= (255.F * vis_brightness);
      col_pixel_recv[ n ].stress_g *= (255.F * vis_brightness);
      col_pixel_recv[ n ].stress_b *= (255.F * vis_brightness);
      
      col_pixel_recv[ n ].density = fmaxf(0.F, (col_pixel_recv[ n ].density - vis_density_threshold_min)) * vis_density_threshold_minmax_inv;
      col_pixel_recv[ n ].stress  = col_pixel_recv[ n ].stress * vis_stress_threshold_max_inv;
    }
  
  if (!write_image) return;
  
  image_file = fopen (image_file_name, "w");
  xdrstdio_create (&xdr_image_file, image_file, XDR_ENCODE);
  
  xdr_int (&xdr_image_file, &screen.pixels_x);
  xdr_int (&xdr_image_file, &screen.pixels_y);
  xdr_int (&xdr_image_file, &col_pixels);
  
  for (n = 0; n < col_pixels; n++)
    {
      xdrWritePixel (&col_pixel_recv[ n ], &xdr_image_file, ColourPalette);
    }
  xdr_destroy (&xdr_image_file);
  fclose (image_file);
}


void visCalculateMouseFlowField (ColPixel *col_pixel_p, LBM *lbm)
{
  vis_mouse_pressure = lbmConvertPressureToPhysicalUnits (col_pixel_p->density * Cs2, lbm);
  
  vis_stess_pressure = lbmConvertStressToPhysicalUnits (col_pixel_p->stress, lbm);
}


void visEnd (void)
{
  rtEnd ();
  
  free(col_pixel_id);
  //free(col_pixel_lock);
  free(col_pixel_recv);
  // free(col_pixel_send);
}

