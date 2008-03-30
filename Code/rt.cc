#include "config.h"


void (*rtAABBvsRay[2][2][2]) (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t_near, float *t_far);

int (*rtTraverseBlocks[2][2][2]) (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]));

int (*rtTraverseVoxels000[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));
int (*rtTraverseVoxels001[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));
int (*rtTraverseVoxels010[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));
int (*rtTraverseVoxels011[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));
int (*rtTraverseVoxels100[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));
int (*rtTraverseVoxels101[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));
int (*rtTraverseVoxels110[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));
int (*rtTraverseVoxels111[2]) (float block_min[], float block_x[], float voxel_flow_field[], float t,
			       void (*ColourPalette) (float value, float col[]));

void (*rtRayTracing[2]) (void (*ColourPalette) (float value, float col[]));


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
#ifndef NOSIMD
  col[3] += dt * palette[3];
#endif
}


int rtTraverseVoxelsVR000 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS000 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsVR001 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[2] += ray_inv[2];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS001 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[2] += ray_inv[2];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsVR010 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[1] += ray_inv[1];
 
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS010 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[1] += ray_inv[1];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsVR011 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[0] -= ray_inv[0];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS011 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[0] -= ray_inv[0];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= block_size2) < 0) return 0;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsVR100 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[0] += ray_inv[0];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS100 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  t_max[0] += ray_inv[0];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsVR101 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[1] -= ray_inv[1];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS101 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[1] -= ray_inv[1];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= block_size) < 0) return 0;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsVR110 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[2] -= ray_inv[2];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS110 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  t_max[2] -= ray_inv[2];
  
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsVR111 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float palette[4];
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  float value;
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      value = voxel_flow_field[ i + j + k ];
      
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[0] - t, palette, ray_col);
		}
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[1] - t, palette, ray_col);
		}
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis_flow_field_cutoff)
		{
		  vis_t_min = t;
		  ColourPalette (value * vis_flow_field_value_max_inv, palette);
		  rtUpdateColour (t_max[2] - t, palette, ray_col);
		}
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

int rtTraverseVoxelsIS111 (float block_min[], float block_x[], float voxel_flow_field[], float t,
			   void (*ColourPalette) (float value, float col[]))
{
  float t_max[VIS_SIMD_SIZE];
  
  int i_vec[VIS_SIMD_SIZE];
  
  int i, j, k;
  
  
 for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (int)block_x[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      i_vec[i] = (i_vec[i] < 0) ? 0 : i_vec[i];
      i_vec[i] = (i_vec[i] > block_size_1) ? block_size_1 : i_vec[i];
    }
  for (i = 0; i < VIS_SIMD_SIZE; i++)
    {
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  i = i_vec[0] * block_size2;
  j = i_vec[1] * block_size;
  k = i_vec[2];
  
  for (;;)
    {
      if ((vis_value = voxel_flow_field[ i + j + k ]) > vis_flow_field_cutoff)
	{
	  vis_t_min = t;
	  return 1;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += block_size2) >= block_size3) return 0;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += block_size) >= block_size2) return 0;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (++k >= block_size) return 0;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}


int rtTraverseBlocks000 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#endif
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels000[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return 0;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels000[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels000[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return 0;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels000[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels000[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

int rtTraverseBlocks001 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#endif
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels001[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  t_max[2] += t_delta[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return 0;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels001[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels001[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return 0;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels001[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels001[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

int rtTraverseBlocks010 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#endif
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels010[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  t_max[1] += t_delta[1];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return 0;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels010[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels010[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return 0;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels010[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels010[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

int rtTraverseBlocks011 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#endif
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels011[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  t_max[0] -= t_delta[0];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i -= cluster_blocks_yz) < 0) return 0;
	      
	      block_min[0] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels011[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels011[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return 0;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels011[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels011[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

int rtTraverseBlocks100 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#endif
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels100[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  t_max[0] += t_delta[0];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return 0;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels100[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels100[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return 0;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels100[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels100[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

int rtTraverseBlocks101 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#endif
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels101[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  t_max[1] -= t_delta[1];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return 0;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels101[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels101[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j -= cluster_blocks_z) < 0) return 0;
	      
	      block_min[1] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels101[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels101[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

int rtTraverseBlocks110 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      i_vec[l] = max(0, min(cluster_blocks_vec[l], (int)(block_size_inv * dx[l])));
      block_min[l] = (float)i_vec[l] * block_size_f - dx[l];
    }
#endif
  i = i_vec[0] * cluster_blocks_yz;
  j = i_vec[1] * cluster_blocks_z;
  k = i_vec[2];
  
  if (block_flow_field[ i+j+k ] != NULL)
    {
      block_x[0] = -block_min[0];
      block_x[1] = -block_min[1];
      block_x[2] = -block_min[2];
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels110[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  t_max[2] -= t_delta[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return 0;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels110[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels110[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return 0;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels110[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (--k < 0) return 0;
	      
	      block_min[2] -= block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels110[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

int rtTraverseBlocks111 (float ray_dx[], float **block_flow_field, void (*ColourPalette) (float value, float col[]))
{
  float block_min[4];
  float t_max[4];
  float block_x[4];
  float t_delta[4];
  float dx[4];
  
  int i_vec[4];
  int i, j, k, l;
  
  
  for (i = 0; i < 3; i++)
    {
      dx[i] = ray_dx[i];
    }
  for (l = 0; l < 4; l++)
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
#ifndef NOSIMD
      block_x[3] = -block_min[3];
#endif
      if (rtTraverseVoxels111[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], 0.F, ColourPalette))
	{
	  return 1;
	}
    }
#ifndef NOSIMD
  for (l = 0; l < 4; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#else
  for (l = 0; l < 3; l++)
    {
      t_max[l] = (block_min[l] + block_size_f) * ray_inv[l];
      t_delta[l] = block_size_f * ray_inv[l];
    }
#endif
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if ((i += cluster_blocks_yz) >= cluster_blocks) return 0;
	      
	      block_min[0] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[0] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[0] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels111[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[0], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels111[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if ((j += cluster_blocks_z) >= cluster_blocks_yz) return 0;
	      
	      block_min[1] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[1] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[1] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels111[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[1], ColourPalette))
		    {
		      return 1;
		    }
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (++k >= cluster_blocks_z) return 0;
	      
	      block_min[2] += block_size_f;
	      
	      if (block_flow_field[ i+j+k ] != NULL)
		{
		  block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
		  block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
		  block_x[2] = t_max[2] * ray_dir[2] - block_min[2];
#ifndef NOSIMD
		  block_x[3] = t_max[2] * ray_dir[3] - block_min[3];
#endif
		  if (rtTraverseVoxels111[vis_mode] (block_min, block_x, block_flow_field[ i+j+k ], t_max[2], ColourPalette))
		    {
		      return 1;
		    }
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
  n = -1;
  
  for (i = 0; i < blocks_x; i++)
    {
      for (j = 0; j < blocks_y; j++)
	{
	  for (k = 0; k < blocks_z; k++)
	    {
	      if (is_block_visited[ ++n ] ||
		  (proc_block_p = &net->proc_block[ n ])->proc_id == NULL)
		{
		  continue;
		}
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
  
  
  cluster_voxel = (float **)malloc(sizeof(float *) * net->my_sites);
  
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
		  
		  cluster_flow_field[ cluster_id ][n] = (float *)malloc(sizeof(float) * sites_in_a_block);
		  
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
				  cluster_flow_field[ cluster_id ][n][m] = -1.F;
				}
			      else
				{
				  cluster_flow_field[ cluster_id ][n][m] = 1.F;
				  cluster_voxel[ my_site_id ] = &cluster_flow_field[ cluster_id ][n][m];
				  
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
  
  rtTraverseVoxels000[0] = rtTraverseVoxelsVR000;
  rtTraverseVoxels001[0] = rtTraverseVoxelsVR001;
  rtTraverseVoxels010[0] = rtTraverseVoxelsVR010;
  rtTraverseVoxels011[0] = rtTraverseVoxelsVR011;
  rtTraverseVoxels100[0] = rtTraverseVoxelsVR100;
  rtTraverseVoxels101[0] = rtTraverseVoxelsVR101;
  rtTraverseVoxels110[0] = rtTraverseVoxelsVR110;
  rtTraverseVoxels111[0] = rtTraverseVoxelsVR111;
  
  rtTraverseVoxels000[1] = rtTraverseVoxelsIS000;
  rtTraverseVoxels001[1] = rtTraverseVoxelsIS001;
  rtTraverseVoxels010[1] = rtTraverseVoxelsIS010;
  rtTraverseVoxels011[1] = rtTraverseVoxelsIS011;
  rtTraverseVoxels100[1] = rtTraverseVoxelsIS100;
  rtTraverseVoxels101[1] = rtTraverseVoxelsIS101;
  rtTraverseVoxels110[1] = rtTraverseVoxelsIS110;
  rtTraverseVoxels111[1] = rtTraverseVoxelsIS111;
  
  rtRayTracing[0] = rtRayTracingVR;
  rtRayTracing[1] = rtRayTracingIS;
  
  rtBuildClusters (net);
}


void rtRayTracingVR (void (*ColourPalette) (float value, float col[]))
{
  // the volume rendering is performed here
  
  float screen_max[4];
  float screen_vtx[4];
  float p0[4], p1[4], p2[4];
  float ray_dx[4];
  float cluster_x[4];
  float dir[4];
  float par1[4], par2[4], par3[4];
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
      
      for (l = 0; l < 4; l++)
	{
	  par3[l] = screen_vtx[l] + subimage_pix[0] * par1[l] + subimage_pix[2] * par2[l];
	}
      for (i = subimage_pix[0]; i <= subimage_pix[1]; i++)
	{
	  for (l = 0; l < 4; l++)
	    {
	      dir[l] = par3[l];
	    }
	  for (j = subimage_pix[2]; j <= subimage_pix[3]; j++)
	    {
	      ray_dir[0] = dir[0];
	      ray_dir[1] = dir[1];
	      ray_dir[2] = dir[2];
#ifndef NOSIMD
	      ray_dir[3] = dir[3];
#endif
	      temp1 = 1.F / sqrtf(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
	      
	      ray_dir[0] *= temp1;
	      ray_dir[1] *= temp1;
	      ray_dir[2] *= temp1;
#ifndef NOSIMD
	      ray_dir[3] *= temp1;
#endif
	      ray_inv[0] = 1.F / ray_dir[0];
	      ray_inv[1] = 1.F / ray_dir[1];
	      ray_inv[2] = 1.F / ray_dir[2];
#ifndef NOSIMD
	      ray_inv[3] = 1.F / ray_dir[3];
#endif
	      ray_sign[0] = ray_dir[0] > 0.F;
	      ray_sign[1] = ray_dir[1] > 0.F;
	      ray_sign[2] = ray_dir[2] > 0.F;
	      
	      dir[0] += par2[0];
	      dir[1] += par2[1];
	      dir[2] += par2[2];
#ifndef NOSIMD
	      dir[3] += par2[3];
#endif
	      // t_near = 0.F;
	      
	      // if (!viewpoint_flag)
		{
		  (*rtAABBvsRay[ray_sign[0]][ray_sign[1]][ray_sign[2]])
		    (&aabb, ray_inv[0], ray_inv[1], ray_inv[2], &t_near, &t_far);
		  
		  if (t_near > t_far) continue;
		  
		  ray_dx[0] = t_near * ray_dir[0] - cluster_x[0];
		  ray_dx[1] = t_near * ray_dir[1] - cluster_x[1];
		  ray_dx[2] = t_near * ray_dir[2] - cluster_x[2];
#ifndef NOSIMD
		  ray_dx[3] = t_near * ray_dir[3] - cluster_x[3];
#endif
		}
	      // else
	      // 	{
	      // 	  ray_dx[0] = cluster_x[0];
	      // 	  ray_dx[1] = cluster_x[1];
	      // 	  ray_dx[2] = cluster_x[2];
	      // 	}
	      ray_col[0] = 0.F;
	      ray_col[1] = 0.F;
	      ray_col[2] = 0.F;
#ifndef NOSIMD
	      ray_col[3] = 0.F;
#endif
	      vis_t_min = 1.e+30F;
	      
	      (*rtTraverseBlocks[ray_sign[0]][ray_sign[1]][ray_sign[2]])
		(ray_dx, block_flow_field, ColourPalette);
	      
	      if (vis_t_min >= 1.e+30F) continue;
	      
	      col_pixel.r = ray_col[ 0 ];
	      col_pixel.g = ray_col[ 1 ];
	      col_pixel.b = ray_col[ 2 ];
	      col_pixel.i = i * (1 << 16) + j;
	      
	      visWritePixel (&col_pixel);
	    }
	  par3[0] += par1[0];
	  par3[1] += par1[1];
	  par3[2] += par1[2];
#ifndef NOSIMD
	  par3[3] += par1[3];
#endif
	}
    }
}


void rtRayTracingIS (void (*ColourPalette) (float value, float col[]))
{
  // the iso-surface is performed
  // here
  
  float screen_max[4];
  float screen_vtx[4];
  float p0[4], p1[4], p2[4];
  float ray_dx[4];
  float cluster_x[4];
  float dir[4];
  float par1[4], par2[4], par3[4];
  float subimage_vtx[4];
  float scale_vec[4];
  float t_near, t_far;
  float **block_flow_field;
  //float temp1;
  
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
      
      for (l = 0; l < 4; l++)
	{
	  par3[l] = screen_vtx[l] + subimage_pix[0] * par1[l] + subimage_pix[2] * par2[l];
	}
      for (i = subimage_pix[0]; i <= subimage_pix[1]; i++)
	{
	  for (l = 0; l < 4; l++)
	    {
	      dir[l] = par3[l];
	    }
	  for (j = subimage_pix[2]; j <= subimage_pix[3]; j++)
	    {
	      ray_dir[0] = dir[0];
	      ray_dir[1] = dir[1];
	      ray_dir[2] = dir[2];
#ifndef NOSIMD
	      ray_dir[3] = dir[3];
#endif
	      ray_inv[0] = 1.F / ray_dir[0];
	      ray_inv[1] = 1.F / ray_dir[1];
	      ray_inv[2] = 1.F / ray_dir[2];
#ifndef NOSIMD
	      ray_inv[3] = 1.F / ray_dir[3];
#endif
	      ray_sign[0] = ray_dir[0] > 0.F;
	      ray_sign[1] = ray_dir[1] > 0.F;
	      ray_sign[2] = ray_dir[2] > 0.F;
	      
	      dir[0] += par2[0];
	      dir[1] += par2[1];
	      dir[2] += par2[2];
#ifndef NOSIMD
	      dir[3] += par2[3];
#endif
	      // t_near = 0.F;
	      
	      // if (!viewpoint_flag)
		{
		  (*rtAABBvsRay[ray_sign[0]][ray_sign[1]][ray_sign[2]])
		    (&aabb, ray_inv[0], ray_inv[1], ray_inv[2], &t_near, &t_far);
		  
		  if (t_near > t_far) continue;
		  
		  ray_dx[0] = t_near * ray_dir[0] - cluster_x[0];
		  ray_dx[1] = t_near * ray_dir[1] - cluster_x[1];
		  ray_dx[2] = t_near * ray_dir[2] - cluster_x[2];
#ifndef NOSIMD
		  ray_dx[3] = t_near * ray_dir[3] - cluster_x[3];
#endif
		}
//	      else
//		{
//		  ray_dx[0] = cluster_x[0];
//		  ray_dx[1] = cluster_x[1];
//		  ray_dx[2] = cluster_x[2];
//#ifndef NOSIMD
//		  ray_dx[3] = cluster_x[3];
//#endif
//		}
	      vis_t_min = 1.e+30F;
	      
	      (*rtTraverseBlocks[ray_sign[0]][ray_sign[1]][ray_sign[2]])
		(ray_dx, block_flow_field, ColourPalette);
	      
	      if (vis_t_min >= 1.e+30F) continue;
	      
	      col_pixel.v = vis_value;
	      col_pixel.t = vis_t_min + t_near;
	      col_pixel.i = i * (1 << 16) + j;
	      
	      visWritePixel (&col_pixel);
	    }
	  par3[0] += par1[0];
	  par3[1] += par1[1];
	  par3[2] += par1[2];
#ifndef NOSIMD
	  par3[3] += par1[3];
#endif
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
  float x1[VIS_SIMD_SIZE], x2[VIS_SIMD_SIZE];
  float temp;
  
  int l;
  
  
  for (l = 0; l < VIS_SIMD_SIZE; l++)
    {
      x1[l] = p1[l] - viewpoint.x[l];
    }
  temp = viewpoint.cos_1 * x1[2] + viewpoint.sin_1 * x1[0];
  
  x2[0] = viewpoint.cos_1 * x1[0] - viewpoint.sin_1 * x1[2];
  x2[1] = viewpoint.cos_2 * x1[1] - viewpoint.sin_2 * temp;
  x2[2] = viewpoint.cos_2 * temp + viewpoint.sin_2 * x1[1];
  
  temp = viewpoint.dist / (-x2[2]);
  
  for (l = 0; l < VIS_SIMD_SIZE; l++)
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
#ifndef NOMPI
  int err;
#endif
  
  
  i = PixelI(col_pixel_p->i);
  j = PixelJ(col_pixel_p->i);
  
  if (*(col_pixel_id_p = &col_pixel_id[ i * screen.pixels_y + j ]) == -1)
    {
      if (col_pixels >= COLOURED_PIXELS_PER_PROC_MAX)
	{
	  printf (" too many coloured pixels per proc\n");
	  printf (" the execution is terminated\n");
#ifndef NOMPI
	  err = MPI_Abort (MPI_COMM_WORLD, 1);
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
      
      if (vis_mode == 0)
	{
	  col_pixel_send[ col_pixels ].r = col_pixel_p->r;
	  col_pixel_send[ col_pixels ].g = col_pixel_p->g;
	  col_pixel_send[ col_pixels ].b = col_pixel_p->b;
	}
      else
	{
	  col_pixel_send[ col_pixels ].v = col_pixel_p->v;
	  col_pixel_send[ col_pixels ].t = col_pixel_p->t;
	}
      col_pixel_send[ col_pixels ].i = col_pixel_p->i;
      ++col_pixels;
    }
  else
    {
      if (vis_mode == 0)
	{
	  col_pixel_send[ *col_pixel_id_p ].r += col_pixel_p->r;
	  col_pixel_send[ *col_pixel_id_p ].g += col_pixel_p->g;
	  col_pixel_send[ *col_pixel_id_p ].b += col_pixel_p->b;
	}
      else
	{
	  if (col_pixel_p->t < col_pixel_send[ *col_pixel_id_p ].t)
	    {
	      col_pixel_send[ *col_pixel_id_p ].v = col_pixel_p->v;
	      col_pixel_send[ *col_pixel_id_p ].t = col_pixel_p->t;
	    }
	}
    }
}


#ifdef STEER
void visReadParameters (char *parameters_file_name, Net *net, Vis *vis, SteerParams *steer)
#else
void visReadParameters (char *parameters_file_name, Net *net, Vis *vis)
#endif
{
  FILE *parameters_file;
  
  float par_to_send[14];
  float ctr_x, ctr_y, ctr_z;
  float longitude, latitude;
  float zoom;
  float density_max, velocity_max, stress_max;
  float dummy;
  
  
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
      
      fscanf (parameters_file, "%e \n", &dummy);
      fscanf (parameters_file, "%e \n", &dummy);
      fscanf (parameters_file, "%e \n", &ctr_x);
      fscanf (parameters_file, "%e \n", &ctr_y);
      fscanf (parameters_file, "%e \n", &ctr_z);
      fscanf (parameters_file, "%e \n", &longitude);
      fscanf (parameters_file, "%e \n", &latitude);
      fscanf (parameters_file, "%e \n", &zoom);
      
      fscanf (parameters_file, "%i \n", &vis_image_freq);
      fscanf (parameters_file, "%i \n", &vis_flow_field_type);
      fscanf (parameters_file, "%i \n", &vis_mode);
      fscanf (parameters_file, "%e \n", &vis_absorption_factor);
      fscanf (parameters_file, "%e \n", &vis_cutoff);
      fscanf (parameters_file, "%e \n", &density_max);
      fscanf (parameters_file, "%e \n", &velocity_max);
      fscanf (parameters_file, "%e \n", &stress_max);

      fclose (parameters_file);
      
      par_to_send[  0 ] = ctr_x;
      par_to_send[  1 ] = ctr_y;
      par_to_send[  2 ] = ctr_z;
      par_to_send[  3 ] = longitude;
      par_to_send[  4 ] = latitude;
      par_to_send[  5 ] = zoom;
      par_to_send[  6 ] = 0.1 + (float)vis_image_freq;
      par_to_send[  7 ] = 0.1 + (float)vis_flow_field_type;
      par_to_send[  8 ] = 0.1 + (float)vis_mode;
      par_to_send[  9 ] = vis_absorption_factor;
      par_to_send[ 10 ] = vis_cutoff;
      par_to_send[ 11 ] = density_max;
      par_to_send[ 12 ] = velocity_max;
      par_to_send[ 13 ] = stress_max;
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 14, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  ctr_x                  =      par_to_send[  0 ];
  ctr_y                  =      par_to_send[  1 ];
  ctr_z                  =      par_to_send[  2 ];
  longitude              =      par_to_send[  3 ];
  latitude               =      par_to_send[  4 ];
  zoom                   =      par_to_send[  5 ];
  vis_image_freq         = (int)par_to_send[  6 ];
  vis_flow_field_type    = (int)par_to_send[  7 ];
  vis_mode               = (int)par_to_send[  8 ];
  vis_absorption_factor  =      par_to_send[  9 ];
  vis_cutoff             =      par_to_send[ 10 ];
  density_max            =      par_to_send[ 11 ];
  velocity_max           =      par_to_send[ 12 ];
  stress_max             =      par_to_send[ 13 ];
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
		 PIXELS_X, PIXELS_Y,
		 ctr_x, ctr_y, ctr_z,
		 5.F * vis->system_size,
		 longitude, latitude,
		 0.5F * (5.F * vis->system_size),
		 zoom);
  
  if (vis_flow_field_type == DENSITY)
    {
      vis_flow_field_value_max_inv = 1.F / density_max;
    }
  else if (vis_flow_field_type == VELOCITY)
    {
      vis_flow_field_value_max_inv = 1.F / velocity_max;
    }
  else
    {
      vis_flow_field_value_max_inv = 1.F / stress_max;
    }
  
#ifdef STEER
  // set up the ReG struct
  
  steer->longitude       = longitude;
  steer->latitude        = latitude;
  steer->zoom            = zoom;
  steer->image_freq      = vis_image_freq;
  steer->flow_field_type = vis_flow_field_type;
  steer->mode            = vis_mode;
  steer->abs_factor      = vis_absorption_factor;
  steer->cutoff          = vis_cutoff;
  steer->max_density     = density_max;
  steer->max_velocity    = velocity_max;
  steer->max_stress      = stress_max;
#endif
}


#ifdef STEER
void visUpdateParameters (Vis *vis, SteerParams *steer)
{
  // update vis params
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
		 PIXELS_X, PIXELS_Y,
		 steer->ctr_x, steer->ctr_y, steer->ctr_z,
		 5.F * vis->system_size,
		 steer->longitude, steer->latitude,
		 0.5F * (5.F * vis->system_size),
		 steer->zoom);
  
  vis_image_freq        = steer->image_freq;
  vis_flow_field_type   = steer->flow_field_type;
  vis_mode              = steer->mode;
  vis_absorption_factor = steer->abs_factor;
  vis_cutoff            = steer->cutoff;
  
  if (vis_flow_field_type == DENSITY)
    {
      vis_flow_field_value_max_inv = 1.F / density_max;
    }
  else if (vis_flow_field_type == VELOCITY)
    {
      vis_flow_field_value_max_inv = 1.F / velocity_max;
    }
  else
    {
      vis_flow_field_value_max_inv = 1.F / stress_max;
    }
}
#endif
 

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
  int col_pixel_count = 5;
  int col_pixel_blocklengths[5] = {1, 1, 1, 1, 1};
  
  MPI_Datatype col_pixel_types[5] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				     MPI_INT,
				     MPI_UB};
  
  MPI_Aint col_pixel_disps[5];
#endif
  
  col_pixels_max = 512 * 512;
  
  // col_pixel_send = (ColPixel *)malloc(sizeof(ColPixel) *  col_pixels_max * max(1, (net_machines - 1)));
  col_pixel_recv = (ColPixel *)malloc(sizeof(ColPixel) * col_pixels_max);
#ifdef RG
  col_pixel_locked = (ColPixel *)malloc(sizeof(ColPixel) * col_pixels_max);
#endif
  
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
  /*
  MPI_Address( &col_pixel_send[0].r, col_pixel_disps + 0 );
  MPI_Address( &col_pixel_send[0].g, col_pixel_disps + 1 );
  MPI_Address( &col_pixel_send[0].b, col_pixel_disps + 2 );
  MPI_Address( &col_pixel_send[0].i, col_pixel_disps + 3 );
  MPI_Address( &col_pixel_send[0]+1, col_pixel_disps + 4 );
  
  int base = col_pixel_disps[0];

  for(int i=0; i<col_pixel_count; i++)
    col_pixel_disps[i] -= base;
  */
  MPI_Type_struct (col_pixel_count, col_pixel_blocklengths, col_pixel_disps, col_pixel_types, &MPI_col_pixel_type);
  MPI_Type_commit (&MPI_col_pixel_type);
#endif
  
  rtInit (net);
}


void visRenderA (void (*ColourPalette) (float value, float col[]), Net *net, Vis *vis)
{
  int pixels_x, pixels_y;
  int i, j;
  int m, n;
  int *col_pixel_id_p;
  int col_pixels_temp;
  int comm_inc, send_id, recv_id;
  int machine_id, master_proc_id;
  
  ColPixel *col_pixel1, *col_pixel2;
  
  
  vis_flow_field_cutoff = vis_cutoff / vis_flow_field_value_max_inv;
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  if (pixels_x * pixels_y > vis_pixels_max)
    {
      vis_pixels_max = pixels_x * pixels_y;
      
      col_pixel_id = (int *)realloc(col_pixel_id, sizeof(int) * vis_pixels_max);
    }
  col_pixels = 0;
  
  (*rtRayTracing[vis_mode]) (ColourPalette);
  
  if (!vis_compositing)
    {
      for (m = 0; m < col_pixels; m++)
	{
	  col_pixel_id[ (PixelI (col_pixel_send[m].i)) * pixels_y + (PixelJ (col_pixel_send[m].i)) ] = -1;
	}
      return;
    }
  
  // here, intra-machine communications are handled through a binary
  // tree pattern and parallel pairwise blocking communications. The
  // master processor of the current machine gets the sub-images of
  // all the processors of that machine. Inter-machine communications,
  // needed if the number of machines is greater than one, take place
  // in the routine visRenderB.
  
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
	      net->err = MPI_Ssend (&col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
	      if (col_pixels > 0)
		{
#ifndef NOMPI
		  net->err = MPI_Ssend (&col_pixel_send, col_pixels, MPI_col_pixel_type,
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
#ifdef RG
			  col_pixel_locked = (ColPixel *)realloc(col_pixel_locked,
								 sizeof(ColPixel) * col_pixels_max);
#endif
			}
		      col_pixel2 = &col_pixel_recv[ *col_pixel_id_p = col_pixels ];
		      
		      memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		      ++col_pixels;
		    }
		  else
		    {
		      col_pixel2 = &col_pixel_recv[ *col_pixel_id_p ];
		      
		      if (vis_mode == 0)
			{
			  col_pixel2->r += col_pixel1->r;
			  col_pixel2->g += col_pixel1->g;
			  col_pixel2->b += col_pixel1->b;
			}
		      else
			{
			  if (col_pixel1->t < col_pixel2->t)
			    {
			      col_pixel2->v = col_pixel1->v;
			      col_pixel2->t = col_pixel1->t;
			    }
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
  
  if (net_machines == 1 || (net->id != 0 && net->id != master_proc_id))
    {
      return;
    }
  
  // inter-machine communications of sub-images begin here
  
  if (net->id != 0)
    {
      recv_id = 0;
#ifndef NOMPI
      net->err = MPI_Ssend (&col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
      if (col_pixels > 0)
	{
	  memcpy (col_pixel_send, col_pixel_recv,
		  col_pixels * sizeof(ColPixel));
#ifndef NOMPI
	  net->err = MPI_Isend (col_pixel_send,
				col_pixels, MPI_col_pixel_type,
				recv_id, 30, MPI_COMM_WORLD,
				&net->req[ 1 ][ net->id * net->procs + recv_id ]);
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
				    &net->req[ 1 ][ (net->id + net->procs) * net->procs + send_id ]);
#endif
	    }
	  send_id += net->procs_per_machine[ m ];
	}
    }
}


void visRenderB (char *image_file_name, void (*ColourPalette) (float value, float col[]), Net *net, Vis *vis)
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
	  for (m = 0; m < col_pixels; m++)
	    {
	      col_pixel_id[ (PixelI (col_pixel_recv[m].i) * pixels_y + PixelJ (col_pixel_recv[m].i)) ] = -1;
	    }
	  return;
	}
      
      if (net->id != 0)
	{
	  recv_id = 0;
#ifndef NOMPI
	  net->err = MPI_Wait (&net->req[ 1 ][ net->id * net->procs + recv_id ], net->status);
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
	      net->err = MPI_Wait (&net->req[ 1 ][ (net->id + net->procs) * net->procs + send_id ], net->status);
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
		      
		      if (vis_mode == 0)
			{
			  col_pixel2->r += col_pixel1->r;
			  col_pixel2->g += col_pixel1->g;
			  col_pixel2->b += col_pixel1->b;
			}
		      else
			{
			  if (col_pixel1->t < col_pixel2->t)
			    {
			      col_pixel2->v = col_pixel1->v;
			      col_pixel2->t = col_pixel1->t;
			    }
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
  
#ifndef BENCH
  
  float factor;
  
  if (net->id != 0) return;
  
  if (vis_mode == 0)
    {
      factor = 255.F * vis_absorption_factor;
    }
  else
    {
      factor = 255.F;
    }
  
#ifdef RG
  
  float col[3];
  
  
  memcpy (col_pixel_locked, col_pixel_recv,
  	  col_pixels * sizeof(ColPixel));
  col_pixels_locked = col_pixels;
  
  if (vis_mode == 0)
    {
      for (n = 0; n < col_pixels; n++)
	{
	  col_pixel_locked[ n ].r *= factor;
	  col_pixel_locked[ n ].g *= factor;
	  col_pixel_locked[ n ].b *= factor;
	}
    }
  else
    {
      for (n = 0; n < col_pixels; n++)
	{
	  ColourPalette (col_pixel_locked[ n ].v * vis_flow_field_value_max_inv, col);
	  
	  col_pixel_locked[ n ].r = col[0] * factor;
	  col_pixel_locked[ n ].g = col[1] * factor;
	  col_pixel_locked[ n ].b = col[2] * factor;
	}
    }
  
#else
    
  FILE *image_file;
  XDR	xdr_image_file;
  
  float voxel_value;
  
  int bits_per_char = sizeof(char) * 8;
  int bits_per_two_chars = 2 * bits_per_char;
  int colour_id, pixel_id;
  int pixel_i, pixel_j;
  
  unsigned char pixel_r, pixel_g, pixel_b;
  
  ColPixel *col_pixel_p;
  
  
  image_file = fopen (image_file_name, "w");
  xdrstdio_create (&xdr_image_file, image_file, XDR_ENCODE);
  
  xdr_int (&xdr_image_file, &vis_mode);
  xdr_int (&xdr_image_file, &vis_flow_field_type);
  
  xdr_int (&xdr_image_file, &screen.pixels_x);
  xdr_int (&xdr_image_file, &screen.pixels_y);
  xdr_int (&xdr_image_file, &col_pixels);
  
  if (vis_mode == 0)
    {
      for (n = 0; n < col_pixels; n++)
	{
	  col_pixel_p = &col_pixel_recv[ n ];
	  
	  pixel_r = (unsigned char)max(0, min(255, (int)(factor * col_pixel_p->r)));
	  pixel_g = (unsigned char)max(0, min(255, (int)(factor * col_pixel_p->g)));
	  pixel_b = (unsigned char)max(0, min(255, (int)(factor * col_pixel_p->b)));
	  
	  pixel_i = PixelI (col_pixel_p->i);
	  pixel_j = PixelJ (col_pixel_p->i);
	  
	  pixel_id = (pixel_i << bits_per_two_chars) + pixel_j;
	  colour_id = (pixel_r << bits_per_two_chars) + (pixel_g << bits_per_char) + pixel_b;
	  
	  xdr_int (&xdr_image_file, &pixel_id);
	  xdr_int (&xdr_image_file, &colour_id);
	}
    }
  else
    {
      for (n = 0; n < col_pixels; n++)
	{
	  col_pixel_p = &col_pixel_recv[ n ];
	  
	  voxel_value = col_pixel_p->v;
	  
	  pixel_i = PixelI (col_pixel_p->i);
	  pixel_j = PixelJ (col_pixel_p->i);
	  
	  pixel_id = (pixel_i << bits_per_two_chars) + pixel_j;
	  
	  xdr_int   (&xdr_image_file, &pixel_id);
	  xdr_float (&xdr_image_file, &voxel_value);
	}
    }
  xdr_destroy (&xdr_image_file);
  fclose (image_file);
  
#endif // RG
#endif // BENCH
}


void visEnd (void)
{
  rtEnd ();
  
  free(col_pixel_id);
#ifdef RG
  free(col_pixel_locked);
#endif
  free(col_pixel_recv);
  // free(col_pixel_send);
}

