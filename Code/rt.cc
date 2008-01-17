#include "config.h"


void (*rtRayAABBIntersection[8]) (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t);

void (*rtTraverse[8]) (float org[],
		       void (*rtAbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		       Cluster *cluster_p, Net *net, Vis *vis);


void rtRayAABBIntersection000 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}

void rtRayAABBIntersection001 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}


void rtRayAABBIntersection010 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}

void rtRayAABBIntersection011 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_1 * inv_x;
  tx1 = aabb->acc_2 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}

void rtRayAABBIntersection100 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}

void rtRayAABBIntersection101 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_3 * inv_y;
  ty1 = aabb->acc_4 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}

void rtRayAABBIntersection110 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_5 * inv_z;
  tz1 = aabb->acc_6 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}

void rtRayAABBIntersection111 (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t)
{
  float tx0, ty0, tz0;
  float tx1, ty1, tz1;
  
  tx0 = aabb->acc_2 * inv_x;
  tx1 = aabb->acc_1 * inv_x;
  ty0 = aabb->acc_4 * inv_y;
  ty1 = aabb->acc_3 * inv_y;
  tz0 = aabb->acc_6 * inv_z;
  tz1 = aabb->acc_5 * inv_z;
  
  if ((*t = fmaxf(tx0, fmaxf(ty0, tz0))) >= fminf(tx1, fminf(ty1, tz1)))
    {
      *t = 1.e+30F;
    }
}


void rtTraverse000 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
	      	{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

void rtTraverse001 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  t_max[2] += ray_inv[2];
  
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

void rtTraverse010 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  t_max[1] += ray_inv[1];
  
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

void rtTraverse011 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  t_max[0] -= ray_inv[0];
  
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i -= block_size2) < 0) return;
	      
	      t = t_max[0];
	      t_max[0] -= ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

void rtTraverse100 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)i_vec[i]) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  t_max[0] += ray_inv[0];
  
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

void rtTraverse101 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  t_max[1] -= ray_inv[1];
  
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		      
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= block_size) < 0) return;
	      
	      t = t_max[1];
	      t_max[1] -= ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

void rtTraverse110 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  t_max[2] -= ray_inv[2];
  
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		      
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max[2];
	      t_max[2] -= ray_inv[2];
	    }
	}
    }
}

void rtTraverse111 (float block_min[], int block_i[], DataBlock *block_p, float t,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float col[]), Vis *vis)
{
  float t_max[VIS_VEC_SIZE];
  float value;
  
  int i, j, k, i_vec[VIS_VEC_SIZE];
  
  unsigned int site_data;
  
  
///#pragma vector always
  for (i = 0; i < VIS_VEC_SIZE; i++)
    {
      i_vec[i] = max(0, min(block_size_1, block_i[i]));
      t_max[i] = (block_min[i] + (float)(i_vec[i] + 1)) * ray_inv[i];
    }
  for (i = 0; i < 2; i++)
    {
      i_vec[i] *= block_size_vec[i];
    }
  i = i_vec[0];
  j = i_vec[1];
  k = i_vec[2];
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[0], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += block_size2) >= block_size3) return;
	      
	      t = t_max[0];
	      t_max[0] += ray_inv[0];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[1], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += block_size) >= block_size2) return;
	      
	      t = t_max[1];
	      t_max[1] += ray_inv[1];
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max[2], ray_col);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= block_size) return;
	      
	      t = t_max[2];
	      t_max[2] += ray_inv[2];
	    }
	}
    }
}

void rtTraverse000 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse000 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (--i_vec[0] < 0) return;
	      ii -= blocks_yz;
	      
	      block_min[0] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse000 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse000 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (--i_vec[1] < 0) return;
	      jj -= blocks_z;
	      
	      block_min[1] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse000 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse000 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverse001 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse001 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  t_max[2] += t_delta[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (--i_vec[0] < 0) return;
	      ii -= blocks_yz;
	      
	      block_min[0] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse001 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse001 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (--i_vec[1] < 0) return;
	      jj -= blocks_z;
	      
	      block_min[1] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse001 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse001 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

void rtTraverse010 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse010 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  t_max[1] += t_delta[1];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (--i_vec[0] < 0) return;
	      ii -= blocks_yz;
	      
	      block_min[0] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse010 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse010 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (++i_vec[1] >= cluster_blocks_vec[1]) return;
	      jj += blocks_z;
	      
	      block_min[1] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse010 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse010 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverse011 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse011 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  t_max[1] += t_delta[1];
  t_max[2] += t_delta[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (--i_vec[0] < 0) return;
	      ii -= blocks_yz;
	      
	      block_min[0] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse011 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] -= t_delta[0];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse011 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (++i_vec[1] >= cluster_blocks_vec[1]) return;
	      jj += blocks_z;
	      
	      block_min[1] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse011 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse011 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

void rtTraverse100 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse100 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  t_max[0] += t_delta[0];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (++i_vec[0] >= cluster_blocks_vec[0]) return;
	      ii += blocks_yz;
	      
	      block_min[0] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse100 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse100 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (--i_vec[1] < 0) return;
	      jj -= blocks_z;
	      
	      block_min[1] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse100 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse100 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverse101 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse101 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  t_max[0] += t_delta[0];
  t_max[2] += t_delta[2];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (++i_vec[0] >= cluster_blocks_vec[0]) return;
	      ii += blocks_yz;
	      
	      block_min[0] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse101 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse101 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (--i_vec[1] < 0) return;
	      jj -= blocks_z;
	      
	      block_min[1] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse101 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] -= t_delta[1];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse101 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}

void rtTraverse110 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse110 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = block_min[l] * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  t_max[0] += t_delta[0];
  t_max[1] += t_delta[1];
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (++i_vec[0] >= cluster_blocks_vec[0]) return;
	      ii += blocks_yz;
	      
	      block_min[0] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse110 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse110 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (++i_vec[1] >= cluster_blocks_vec[1]) return;
	      jj += blocks_z;
	      
	      block_min[1] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse110 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (--i_vec[2] < 0) return;
	      --kk;
	      
	      block_min[2] -= block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse110 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] -= t_delta[2];
	    }
	}
    }
}

void rtTraverse111 (float org[],
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2, float col[]),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float block_min[VIS_VEC_SIZE];
  float t_max[VIS_VEC_SIZE];
  float t_delta[VIS_VEC_SIZE];
  float dx[VIS_VEC_SIZE];
  
  int block_i[VIS_VEC_SIZE];
  int i_vec[VIS_VEC_SIZE], ii_vec[3];
  int cluster_blocks_vec[3];
  int l;
  int ii, jj, kk;

  DataBlock *block_p;
  
  
  cluster_blocks_vec[0] = cluster_p->blocks_x;
  cluster_blocks_vec[1] = cluster_p->blocks_y;
  cluster_blocks_vec[2] = cluster_p->blocks_z;
  
///#pragma vector always
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      dx[l] = org[l] - cluster_p->x[l];
      i_vec[l] = max(0, min(cluster_blocks_vec[l] - 1, (int)(block_size_inv * dx[l])));
    }
  for (l = 0; l < 3; l++)
    {
      ii_vec[l] = (i_vec[l] + cluster_p->block_min[l]) * blocks_vec[l];
      block_min[l] = (float)(i_vec[l] * block_size) - dx[l];
    }
  ii = ii_vec[0];
  jj = ii_vec[1];
  kk = ii_vec[2];
  
  if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      for (l = 0; l < 3; l++)
	{
	  block_i[l] = (int)(-block_min[l]);
	}
      rtTraverse111 (block_min, block_i, block_p, 0.F, AbsorptionCoefficients, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      t_max[l] = (block_min[l] + block_size) * ray_inv[l];
      t_delta[l] = block_size * ray_inv[l];
    }
  
  for (;;)
    {
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      if (++i_vec[0] >= cluster_blocks_vec[0]) return;
	      ii += blocks_yz;
	      
	      block_min[0] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[0] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse111 (block_min, block_i, block_p, t_max[0], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse111 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      if (++i_vec[1] >= cluster_blocks_vec[1]) return;
	      jj += blocks_z;
	      
	      block_min[1] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[1] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse111 (block_min, block_i, block_p, t_max[1], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      if (++i_vec[2] >= cluster_blocks_vec[2]) return;
	      ++kk;
	      
	      block_min[2] += block_size;
	      
	      if ((block_p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  for (l = 0; l < 3; l++)
		    {
		      block_i[l] = (int)(t_max[2] * ray_dir[l] - block_min[l]);
		    }
		  rtTraverse111 (block_min, block_i, block_p, t_max[2], AbsorptionCoefficients, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max[2] += t_delta[2];
	    }
	}
    }
}


void rtInit (Net *net, Vis *vis)
{
  rtRayAABBIntersection[0] = rtRayAABBIntersection000;
  rtRayAABBIntersection[1] = rtRayAABBIntersection001;
  rtRayAABBIntersection[2] = rtRayAABBIntersection010;
  rtRayAABBIntersection[3] = rtRayAABBIntersection011;
  rtRayAABBIntersection[4] = rtRayAABBIntersection100;
  rtRayAABBIntersection[5] = rtRayAABBIntersection101;
  rtRayAABBIntersection[6] = rtRayAABBIntersection110;
  rtRayAABBIntersection[7] = rtRayAABBIntersection111;
  
  rtTraverse[0] = rtTraverse000;
  rtTraverse[1] = rtTraverse001;
  rtTraverse[2] = rtTraverse010;
  rtTraverse[3] = rtTraverse011;
  rtTraverse[4] = rtTraverse100;
  rtTraverse[5] = rtTraverse101;
  rtTraverse[6] = rtTraverse110;
  rtTraverse[7] = rtTraverse111;
}


void rtRayTracing (void (*AbsorptionCoefficients) (float flow_field_data, float t1, float t2, float col[]),
		   Net *net, Vis *vis)
{
  // here, the ray tracing is performed and the intra-machine communications take place
  
  float screen_max[4];
  float screen_vtx[VIS_VEC_SIZE];
  float p0[VIS_VEC_SIZE], p1[VIS_VEC_SIZE], p2[VIS_VEC_SIZE];
  float dir[VIS_VEC_SIZE];
  float cluster_blocks_vec[VIS_VEC_SIZE];
  float t;
  float temp1, temp2;
  
  float par1[VIS_VEC_SIZE], par2[VIS_VEC_SIZE], par3[VIS_VEC_SIZE];
  float subimage_vtx[4];
  float v[2][VIS_VEC_SIZE];
  float scale_vec[4];
  
  int pixels_x, pixels_y;
  int i, j, k, l;
  int cluster_id;
  int subimage_pix[4];
  int viewpoint_flag;
  int ray_dir_code;
  
  AABB aabb;
  
  Cluster *cluster_p;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  screen_max[0] = screen.max_x;
  screen_max[1] = screen.max_x;
  screen_max[2] = screen.max_y;
  screen_max[3] = screen.max_y;
  
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      p0[l] = viewpoint.x[l];
      screen_vtx[l] = EPSILON + screen.ctr[l] - screen.dir1[l] - screen.dir2[l] - p0[l];
    }
  temp1 = (2.F / (float)pixels_x);
  temp2 = (2.F / (float)pixels_y);
  
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      par1[l] = screen.dir1[l] * temp1;
      par2[l] = screen.dir2[l] * temp2;
    }
  scale_vec[0] = scale_vec[1] = 1.F / screen.par_x;
  scale_vec[2] = scale_vec[3] = 1.F / screen.par_y;

  for (cluster_id = 0; cluster_id < vis->clusters; cluster_id++)
    {
      cluster_p = &vis->cluster[ cluster_id ];
      
      // the image-based projection of the cluster bounding box is calculated here
      
      cluster_blocks_vec[0] = (float)cluster_p->blocks_x;
      cluster_blocks_vec[1] = (float)cluster_p->blocks_y;
      cluster_blocks_vec[2] = (float)cluster_p->blocks_z;
      
      for (l = 0; l < VIS_VEC_SIZE; l++)
	{
	  v[0][l] = cluster_p->x[l];
	  v[1][l] = cluster_p->x[l] + cluster_blocks_vec[l] * (float)block_size;
	}
      subimage_vtx[0] =  1.e+30F;
      subimage_vtx[1] = -1.e+30F;
      subimage_vtx[2] =  1.e+30F;
      subimage_vtx[3] = -1.e+30F;
      
      for (i = 0; i < 2; i++)
	{
	  p1[0] = v[i][0];
	  
	  for (j = 0; j < 2; j++)
	    {
	      p1[1] = v[j][1];
	      
	      for (k = 0; k < 2; k++)
		{
		  p1[2] = v[k][2];
		  
		  visProject (p1, p2);
		  
		  subimage_vtx[0] = fminf(subimage_vtx[0], p2[0]);
		  subimage_vtx[1] = fmaxf(subimage_vtx[1], p2[0]);
		  subimage_vtx[2] = fminf(subimage_vtx[2], p2[1]);
		  subimage_vtx[3] = fmaxf(subimage_vtx[3], p2[1]);
		}
	    }
	}
      for (l = 0; l < 4; l++)
	{
	  subimage_pix[l] = (int)(scale_vec[l] * (subimage_vtx[l] + screen_max[l]));
	}
      if (subimage_pix[0] >= pixels_x || subimage_pix[1] < 0 ||
	  subimage_pix[2] >= pixels_y || subimage_pix[3] < 0)
	{
	  continue;
	}
      subimage_pix[0] = max(subimage_pix[0], 0);
      subimage_pix[1] = min(subimage_pix[1], pixels_x - 1);
      subimage_pix[2] = max(subimage_pix[2], 0);
      subimage_pix[3] = min(subimage_pix[3], pixels_y - 1);
      
      if (p0[0] >= v[0][0] && p0[1] >= v[0][1] && p0[2] >= v[0][2] &&
	  p0[0] <= v[1][0] && p0[1] <= v[1][1] && p0[2] <= v[1][2])
	{
	  viewpoint_flag = 1;
	}
      else
	{
	  viewpoint_flag = 0;
	}
      aabb.acc_1 = v[1][0] - p0[0];
      aabb.acc_2 = v[0][0] - p0[0];
      aabb.acc_3 = v[1][1] - p0[1];
      aabb.acc_4 = v[0][1] - p0[1];
      aabb.acc_5 = v[1][2] - p0[2];
      aabb.acc_6 = v[0][2] - p0[2];
      
      temp1 = (float)subimage_pix[0] + 0.5F;
      temp2 = (float)subimage_pix[2] + 0.5F;
      
      for (l = 0; l < VIS_VEC_SIZE; l++)
	{
	  par3[l] = screen_vtx[l] + temp1 * par1[l] + temp2 * par2[l];
	}
      
      for (i = subimage_pix[0]; i <= subimage_pix[1]; i++)
	{
	  for (l = 0; l < VIS_VEC_SIZE; l++)
	    {
	      dir[l] = par3[l];
	    }
	  for (j = subimage_pix[2]; j <= subimage_pix[3]; j++)
	    {
	      for (l = 0; l < VIS_VEC_SIZE; l++)
		{
		  p1[l] = p0[l];
		  ray_dir[l] = dir[l];
		}
	      temp1 = 1.F / sqrtf(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
	      
	      for (l = 0; l < VIS_VEC_SIZE; l++)
		{	      
		  ray_dir[l] *= temp1;
		  ray_inv[l] = 1.F / ray_dir[l];
		}
	      ray_dir_code = ray_dir[2] > 0.F;
	      ray_dir_code |= (ray_dir[1] > 0.F) << 1;
	      ray_dir_code |= (ray_dir[0] > 0.F) << 2;
	      
	      for (l = 0; l < VIS_VEC_SIZE; l++)
		{
		  dir[l] += par2[l];
		}
	      t = 0.F;
	      
	      if (!viewpoint_flag)
		{
		  (*rtRayAABBIntersection[ ray_dir_code ]) (&aabb, ray_inv[0], ray_inv[1], ray_inv[2], &t);
		  
		  if (t >= 1.e+30F) continue;
		  
		  for (l = 0; l < VIS_VEC_SIZE; l++)
		    {
		      p1[l] += t * ray_dir[l];
		    }
		}
	      
	      for (l = 0; l < VIS_VEC_SIZE; l++)
		{
		  ray_col[l] = 0.F;
		}
	      vis->t_min = 1.e+30F;
	      
	      (*rtTraverse[ ray_dir_code ]) (p1, AbsorptionCoefficients, cluster_p, net, vis);
	      
	      if (vis->t_min >= 1.e+30F) continue;
	      
	      visWritePixel (ray_col, vis->t_min + t, i, j, vis);
	      
	    }
	  for (l = 0; l < VIS_VEC_SIZE; l++)
	    {
	      par3[l] += par1[l];
	    }
	}
    }
}


void rtEnd (Vis *vis)
{
  free(vis->cluster);
}


void slInit (Net *net, Vis *vis)
{
  int boundary_sites_count;
  int site_vec[3];
  int i, j, k;
  int l, m, n;
  
  unsigned int site_id, site_type;
  
  DataBlock *map_block_p;
  
  
  vis->seeds_max = 10000;
  vis->seed = (float *)malloc(sizeof(float) * VIS_VEC_SIZE * vis->seeds_max);
  
  vis->seeds = 0;
  
  boundary_sites_count = 0;
  
  n = -1;
  
  for (i = 0; i < sites_x; i += block_size)
    {
      for (j = 0; j < sites_y; j += block_size)
	{
	  for (k = 0; k < sites_z; k += block_size)
	    {
	      if (net->proc_id[ ++n ] != net->id)
		{
		  continue;
		}
	      map_block_p = &net->map_block[ n ];
	      
	      m = -1;
	      
	      for (site_vec[0] = i; site_vec[0] < i + block_size; site_vec[0]++)
		{
		  for (site_vec[1] = j; site_vec[1] < j + block_size; site_vec[1]++)
		    {
		      for (site_vec[2] = k; site_vec[2] < k + block_size; site_vec[2]++)
			{
			  site_id = map_block_p->site_data[ ++m ];
			  
			  if (site_id & (1U << 31U))
			    {
			      continue;
			    }
			  site_type = net->site_data[ site_id ] & SITE_TYPE_MASK;
			  
			  if (site_type != INLET_TYPE &&
			      site_type != OUTLET_TYPE)
			    {
			      continue;
			    }
			  if ((++boundary_sites_count)%10 != 1)
			    {
			      continue;
			    }
			  if (vis->seeds == vis->seeds_max)
			    {
			      vis->seeds_max <<= 1;
			      vis->seed = (float *)realloc(vis->seed,
							   sizeof(float) * 3 * vis->seeds_max);
			    }
			  for (l = 0; l < VIS_VEC_SIZE; l++)
			    {
			      vis->seed[ VIS_VEC_SIZE*vis->seeds+l ] = (float)site_vec[l] + (0.5F - vis->half_dim[l]);
			    }
			  ++vis->seeds;
			}
		    }
		}
	    }
	}
    }
  vis->streamlines_max = 10000;
  vis->streamline = (float *)malloc(sizeof(float) * 3 * vis->streamlines_max);
}


void slStreamlines (void ColourPalette (float vel_m, float col[]),
		    Net *net, Vis *vis)
{
  double density, v[VIS_VEC_SIZE];
  
  float scale_x, scale_y;
  float vel_m;
  float col[3];
  float x[VIS_VEC_SIZE];
  float z, z1, z2, dz;
  float z_old, z_new;
  float d, incE, incNE;
  
  int pixels_x, pixels_y;
  int cycles;
  int is_inside_my_subdomain;
  int proc_id, neigh_proc_index;
  int site_vec[VIS_VEC_SIZE];
  int b_vec[3];
  int s_vec[3];
  int block_id;
  int i_old, j_old;
  int i_new, j_new;
  int i1, j1, i2, j2;
  int di, dj;
  int i, j;
  int l, m, n;
#ifndef NOMPI
  int streamlines_max;
#endif
  unsigned int site_id;
  
  DataBlock *map_block_p;
#ifndef NOMPI
  NeighProc *neigh_proc_p;
#endif
  
  vis->streamlines = vis->seeds;
  memcpy (vis->streamline, vis->seed, vis->seeds * VIS_VEC_SIZE * sizeof(float));
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  scale_x = 1.F / screen.par_x;
  scale_y = 1.F / screen.par_y;
  
  for (cycles = 0; cycles < net->procs; cycles++)
    {
      for (m = 0; m < net->neigh_procs; m++)
	{
	  streamlines_to_send[ m ] = 0;
	}
      for (n = vis->streamlines - 1; n >= 0; n--)
	{
///#pragma vector always
	  for (l = 0; l < VIS_VEC_SIZE; l++)
	    {
	      site_vec[l] = (int)(vis->streamline[ VIS_VEC_SIZE*n+l ] + vis->half_dim[l]);
	    }
	  if (site_vec[0] < 0 || site_vec[0] >= sites_x ||
	      site_vec[1] < 0 || site_vec[1] >= sites_y ||
	      site_vec[2] < 0 || site_vec[2] >= sites_z)
	    {
	      --vis->streamlines;
	      continue;
	    }
	  for (l = 0; l < 3; l++)
	    {
	      b_vec[l] = site_vec[l] / block_size;
	    }
	  block_id = (b_vec[0] * blocks_y + b_vec[1]) * blocks_z + b_vec[2];
	  
	  proc_id = net->proc_id[ block_id ];
	  
	  if (proc_id != net->id)
	    {
	      --vis->streamlines;
	      continue;
	    }
	  for (l = 0; l < 3; l++)
	    {
	      s_vec[l] = site_vec[l] - (b_vec[l] * block_size);
	    }
	  site_id = (((s_vec[0] * block_size) + s_vec[1]) * block_size) + s_vec[2];
	  site_id = net->map_block[ block_id ].site_data[ site_id ];
	  
	  if (site_id & (1U << 31U))
	    {
	      --vis->streamlines;
	      continue;
	    }
	  visProject (&vis->streamline[ VIS_VEC_SIZE*n ], x);
	  
	  i_old = (int)(scale_x * (x[0] + screen.max_x));
	  j_old = (int)(scale_y * (x[1] + screen.max_y));
	  z_old = x[2];
	  
	  is_inside_my_subdomain = 1;
	  
	  while (is_inside_my_subdomain)
	    {
	      is_inside_my_subdomain = 0;
	      
	      lbmDensityAndVelocity (&f_old[ site_id*15 ], &density, &v[0], &v[1], &v[2]);
	      
	      vel_m = (float)sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	      
	      if (vel_m < 1.e-6)
		{
		  --vis->streamlines;
		  continue;
		}
	      ColourPalette (vel_m * vis->flow_field_value_max_inv[ VELOCITY ], col);
	      
	      vel_m = 1.F / vel_m;
///#pragma vector always
	      for (l = 0; l < VIS_VEC_SIZE; l++)
		{
		  vis->streamline[ VIS_VEC_SIZE*n+l ] += (float)v[l] * vel_m;
		}
	      visProject (&vis->streamline[ VIS_VEC_SIZE*n ], x);
	      
	      z = x[2];
	      
	      i_new = (int)(scale_x * (x[0] + screen.max_x));
	      j_new = (int)(scale_y * (x[1] + screen.max_y));
	      z_new = z;
	      
	      i1 = i_old;
	      j1 = j_old;
	      z1 = z_old;
	      i2 = i_new;
	      j2 = j_new;
	      z2 = z_new;
	      
	      if (i2 < i1)
		{
		  i = i1;
		  j = j1;
		  z = z1;
		  i1 = i2;
		  j1 = j2;
		  z1 = z2;
		  i2 = i;
		  j2 = j;
		  z2 = z;
		}
	      di = i2 - i1;
	      if (j1 < j2)
		{
		  m = 1;
		  dj = j2 - j1;
		}
	      else
		{
		  m = -1;
		  dj = j1 - j2;
		}
	      i = i1;
	      j = j1;
	      z = z1;
	      
	      if (di > dj)
		{
		  incE = dj;
		  d = dj - di;
		  incNE = d;
		  dz = (z2 - z1) / (float)max(1, dj);
		  
		  while (i != i2)
		    {
		      if (!(i < 0 || i >= pixels_x ||
			    j < 0 || j >= pixels_y))
			{
			  visWritePixel (col, z, i, j, vis);
			}
		      if (d < 0)
			{
			  d += incE;
			}
		      else
			{
			  d += incNE;
			  j += m;
			  z += dz;
			}
		      ++i;
		    }
		}
	      else
		{
		  incE = di;
		  d = di - dj;
		  incNE = d;
		  dz = (z2 - z1) / (float)max(1, di);
		  
		  while (i != i2)
		    {
		      if (!(i < 0 || i >= pixels_x ||
			    j < 0 || j >= pixels_y))
			{
			  visWritePixel (col, z, i, j, vis);
			}
		      if (d < 0)
			{
			  d += incE;
			}
		      else
			{
			  d += incNE;
			  ++i;
			  z += dz;
			}
		      j += m;
		    }
		}
///#pragma vector always
	      for (l = 0; l < VIS_VEC_SIZE; l++)
		{
		  site_vec[l] = (int)(vis->streamline[ VIS_VEC_SIZE*n+l ] + vis->half_dim[l]);
		}
	      if (site_vec[0] < 0 || site_vec[0] >= sites_x ||
		  site_vec[1] < 0 || site_vec[1] >= sites_y ||
		  site_vec[2] < 0 || site_vec[2] >= sites_z)
		{
		  --vis->streamlines;
		  continue;
		}
	      for (l = 0; l < 3; l++)
		{
		  b_vec[l] = site_vec[l] / block_size;
		}
	      block_id = (b_vec[0] * blocks_y + b_vec[1]) * blocks_z + b_vec[2];
	      
	      proc_id = net->proc_id[ block_id ];
	      
	      if (proc_id != net->id && proc_id != (1 << 14))
		{
		  neigh_proc_index = net->from_proc_id_to_neigh_proc_index[ proc_id ];
		  
		  if (streamlines_to_send[ neigh_proc_index ] < STREAMLINES_MAX)
		    {
		      memcpy (&streamline_to_send[ neigh_proc_index ][ VIS_VEC_SIZE*streamlines_to_send[neigh_proc_index] ],
		  	      &vis->streamline[ VIS_VEC_SIZE*n ], VIS_VEC_SIZE * sizeof(float));
		      ++streamlines_to_send[ neigh_proc_index ];
		    }
		  --vis->streamlines;
		  continue;
		}
	      map_block_p = &net->map_block[ block_id ];
	      
	      if (map_block_p->site_data == NULL)
		{
		  --vis->streamlines;
		  continue;
		}
	      for (l = 0; l < 3; l++)
		{
		  s_vec[l] = site_vec[l] - (b_vec[l] * block_size);
		}
	      site_id = (((s_vec[0] * block_size) + s_vec[1]) * block_size) + s_vec[2];
	      site_id = net->map_block[ block_id ].site_data[ site_id ];
	      
	      if (site_id & (1U << 31U))
		{
		  --vis->streamlines;
		  continue;
		}
	      i_old = i_new;
	      j_old = j_new;
	      z_old = z_new;
	      
	      is_inside_my_subdomain = 1;
	    }
	}
#ifndef NOMPI
      for (m = 0; m < net->neigh_procs; m++)
      	{
      	  neigh_proc_p = &net->neigh_proc[ m ];
      	  
      	  net->err = MPI_Send (&streamlines_to_send[m], 1, MPI_INT,
				neigh_proc_p->id, 30, MPI_COMM_WORLD);
      	  net->err = MPI_Recv (&streamlines_to_recv[m], 1, MPI_INT,
			       neigh_proc_p->id, 30, MPI_COMM_WORLD, net->status);
      	  
	  if (streamlines_to_send[ m ] > 0)
      	    {
      	      net->err = MPI_Send (&streamline_to_send[m][0], VIS_VEC_SIZE * streamlines_to_send[ m ],
				    MPI_FLOAT, neigh_proc_p->id, 30, MPI_COMM_WORLD);
      	    }
	  if (streamlines_to_recv[ m ] > 0)
      	    {
      	      net->err = MPI_Recv (&streamline_to_recv[m][0], VIS_VEC_SIZE * streamlines_to_recv[ m ],
				   MPI_FLOAT, neigh_proc_p->id, 30, MPI_COMM_WORLD, net->status);
	      
              streamlines_max = min(STREAMLINES_MAX - vis->streamlines,
      	  			    streamlines_to_recv[ m ]);
      	      
      	      memcpy (&vis->streamline[ vis->streamlines ], &streamline_to_recv[ m ][ 0 ],
      	  	      streamlines_max * VIS_VEC_SIZE * sizeof(float));
      	      vis->streamlines += streamlines_max;
      	    }
	}
#endif
    }
}


void slEnd (Vis *vis)
{
  free(vis->streamline);
  vis->streamline = NULL;
  
  free(vis->seed);
  vis->seed = NULL;
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
  float x1[VIS_VEC_SIZE], x2[VIS_VEC_SIZE];
  float temp;
  
  int l;
  
  
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      x1[l] = p1[l] - viewpoint.x[l];
    }
  temp = viewpoint.cos_1 * x1[2] + viewpoint.sin_1 * x1[0];
  
  x2[0] = viewpoint.cos_1 * x1[0] - viewpoint.sin_1 * x1[2];
  x2[1] = viewpoint.cos_2 * x1[1] - viewpoint.sin_2 * temp;
  x2[2] = viewpoint.cos_2 * temp + viewpoint.sin_2 * x1[1];
  
  p2[2] = -x2[2];
  
  temp = screen.dist / p2[2];
  
  for (l = 0; l < VIS_VEC_SIZE; l++)
    {
      p2[l] = temp * x2[l];
    }
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
  
  temp = dist / rad;
  
  screen.ctr[0] = viewpoint.x[0] + temp * (ctr_x - viewpoint.x[0]);
  screen.ctr[1] = viewpoint.x[1] + temp * (ctr_y - viewpoint.x[1]);
  screen.ctr[2] = viewpoint.x[2] + temp * (ctr_z - viewpoint.x[2]);
  
  screen.zoom = zoom;
  screen.dist = dist;
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     screen.max_x, 0.0F, 0.0F,
	     &screen.dir1[0], &screen.dir1[1], &screen.dir1[2]);
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     0.0F, screen.max_y, 0.0F,
	     &screen.dir2[0], &screen.dir2[1], &screen.dir2[2]);
  
  screen.par_x = (2.F * screen.max_x) / (float)pixels_x;
  screen.par_y = (2.F * screen.max_y) / (float)pixels_y;
}


void visWritePixel (float col[], float t, int i, int j, Vis *vis)
{
  int *col_pixel_id_p;
#ifndef NOMPI
  int err;
#endif
  
  ColPixel *col_pixel_p;
  
  
  if (*(col_pixel_id_p = &vis->col_pixel_id[ i * screen.pixels_y + j ]) == -1)
    {
      if (vis->col_pixels >= COLOURED_PIXELS_PER_PROC_MAX)
	{
	  printf (" too many coloured pixels per proc\n");
	  printf (" the execution is terminated\n");
#ifndef NOMPI
	  err = MPI_Abort (MPI_COMM_WORLD, 1);
#else
	  exit(1);
#endif
	}
      if (vis->col_pixels == vis->col_pixels_max)
	{
	  vis->col_pixels_max <<= 1;
	  vis->col_pixel_recv = (ColPixel *)realloc(vis->col_pixel_recv,
						    sizeof(ColPixel) * vis->col_pixels_max);
	}
      *col_pixel_id_p = vis->col_pixels;
      
      col_pixel_p = &vis->col_pixel_send[ *col_pixel_id_p ];
      col_pixel_p->r = 0.F;
      col_pixel_p->g = 0.F;
      col_pixel_p->b = 0.F;
      col_pixel_p->t = 1.e+30F;
      col_pixel_p->i = (short int)i;
      col_pixel_p->j = (short int)j;
      ++vis->col_pixels;
    }
  else
    {
      col_pixel_p = &vis->col_pixel_send[ *col_pixel_id_p ];
    }
  if (vis->mode == 0)
    {
      col_pixel_p->r += col[0];
      col_pixel_p->g += col[1];
      col_pixel_p->b += col[2];
    }
  else if (vis->mode == 1 || vis->mode == 2)
    {
      if (t < col_pixel_p->t)
	{
	  col_pixel_p->r = col[0];
	  col_pixel_p->g = col[1];
	  col_pixel_p->b = col[2];
	  col_pixel_p->t = t;
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
      
      fscanf (parameters_file, "%i \n", &vis->image_freq);
      fscanf (parameters_file, "%i \n", &vis->flow_field_type);
      fscanf (parameters_file, "%i \n", &vis->mode);
      fscanf (parameters_file, "%e \n", &vis->absorption_factor);
      fscanf (parameters_file, "%e \n", &vis->cutoff);
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
      par_to_send[  6 ] = 0.1 + (float)vis->image_freq;
      par_to_send[  7 ] = 0.1 + (float)vis->flow_field_type;
      par_to_send[  8 ] = 0.1 + (float)vis->mode;
      par_to_send[  9 ] = vis->absorption_factor;
      par_to_send[ 10 ] = vis->cutoff;
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
  vis->image_freq        = (int)par_to_send[  6 ];
  vis->flow_field_type   = (int)par_to_send[  7 ];
  vis->mode              = (int)par_to_send[  8 ];
  vis->absorption_factor =      par_to_send[  9 ];
  vis->cutoff            =      par_to_send[ 10 ];
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
  
  vis->flow_field_value_max_inv[ DENSITY  ] = 1.F / density_max;
  vis->flow_field_value_max_inv[ VELOCITY ] = 1.F / velocity_max;
  vis->flow_field_value_max_inv[ STRESS   ] = 1.F / stress_max;
  
#ifdef STEER
  // set up the ReG struct
  
  steer->longitude       = longitude;
  steer->latitude        = latitude;
  steer->zoom            = zoom;
  steer->image_freq      = vis->image_freq;
  steer->flow_field_type = vis->flow_field_type;
  steer->mode            = vis->mode;
  steer->abs_factor      = vis->absorption_factor;
  steer->cutoff          = vis->cutoff;
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
  
  vis->image_freq        = steer->image_freq;
  vis->flow_field_type   = steer->flow_field_type;
  vis->mode              = steer->mode;
  vis->absorption_factor = steer->abs_factor;
  vis->cutoff            = steer->cutoff;
  
  vis->flow_field_value_max_inv[ DENSITY  ] = 1.F / steer->max_density;
  vis->flow_field_value_max_inv[ VELOCITY ] = 1.F / steer->max_velocity;
  vis->flow_field_value_max_inv[ STRESS   ] = 1.F / steer->max_stress;
}
#endif
 

void visInit (char *image_file_name, Net *net, Vis *vis)
{
  blocks_yz = blocks_y * blocks_z;
  
  vis->half_dim[0] = 0.5F * (float)sites_x;
  vis->half_dim[1] = 0.5F * (float)sites_y;
  vis->half_dim[2] = 0.5F * (float)sites_z;
  
  vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1], vis->half_dim[2]));
  
  block_size2 = block_size * block_size;
  block_size3 = block_size * block_size2;
  block_size_1 = block_size - 1;
  
  block_size_inv = 1.F / (float)block_size;
  
  blocks_vec[0] = blocks_y * blocks_z;
  blocks_vec[1] = blocks_z;
  blocks_vec[2] = 1;
  blocks_vec[3] = 0;
  
  block_size_vec[0] = block_size * block_size;
  block_size_vec[1] = block_size;
  block_size_vec[2] = 1;
  block_size_vec[3] = 0;
  
  
#ifndef NOMPI
  int col_pixel_count = 7;
  int col_pixel_blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
  
  MPI_Datatype col_pixel_types[7] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				     MPI_SHORT, MPI_SHORT,
				     MPI_UB};
  
  MPI_Aint col_pixel_disps[7];
#endif
  
  
  vis->image_file_name = image_file_name;
  
  vis->col_pixels_max = 512 * 512;
  
  /// vis->col_pixel_send = (ColPixel *)malloc(sizeof(ColPixel) * vis->col_pixels_max);
  
  vis->col_pixel_recv = (ColPixel *)malloc(sizeof(ColPixel) * vis->col_pixels_max);
  vis->col_pixel_locked = (ColPixel *)malloc(sizeof(ColPixel) * vis->col_pixels_max);
  
  vis->pixels_max = IMAGE_SIZE;
  vis->col_pixel_id = (int *)malloc(sizeof(int) * vis->pixels_max);
  
  
#ifndef NOMPI
  // create the derived datatype for the MPI communications
  /*
  col_pixel_disps[0] = 0;
  
  for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(float) * col_pixel_blocklengths[i - 1]);
	}
      else if (col_pixel_types[i - 1] == MPI_SHORT)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(short int) * col_pixel_blocklengths[i - 1]);
	}
    }
  */
  MPI_Address( &vis->col_pixel_send[0].r, col_pixel_disps + 0 ); 
  MPI_Address( &vis->col_pixel_send[0].g, col_pixel_disps + 1 ); 
  MPI_Address( &vis->col_pixel_send[0].b, col_pixel_disps + 2 ); 
  MPI_Address( &vis->col_pixel_send[0].t, col_pixel_disps + 3 ); 
  MPI_Address( &vis->col_pixel_send[0].i, col_pixel_disps + 4 );
  MPI_Address( &vis->col_pixel_send[0].j, col_pixel_disps + 5 );
  MPI_Address( &vis->col_pixel_send[0]+1, col_pixel_disps + 6 );
  
  int base = col_pixel_disps[0];

  for(int i=0; i<col_pixel_count; i++)
    col_pixel_disps[i] -= base;
  
  MPI_Type_struct (col_pixel_count, col_pixel_blocklengths, col_pixel_disps, col_pixel_types, &MPI_col_pixel_type);
  MPI_Type_commit (&MPI_col_pixel_type);
#endif
  
  rtInit (net, vis);
  
  slInit (net, vis);
}


void visRenderA (void (*rtAbsorptionCoefficients) (float flow_field_data, float t1, float t2, float col[]),
		 void slColourPalette (float vel_m, float col[]),
		 Net *net, Vis *vis)
{
  int pixels_x, pixels_y;
  int i, j;
  int m, n;
  int *col_pixel_id_p;
  int col_pixels;
  int comm_inc, send_id, recv_id;
  int machine_id, master_proc_id;
  
  ColPixel *col_pixel1, *col_pixel2;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  if (pixels_x * pixels_y > vis->pixels_max)
    {
      vis->pixels_max = pixels_x * pixels_y;
      
      vis->col_pixel_id = (int *)realloc(vis->col_pixel_id, sizeof(int) * vis->pixels_max);
    }
  for (i = 0; i < pixels_x * pixels_y; i++)
    {
      vis->col_pixel_id[ i ] = -1;
    }
  vis->col_pixels = 0;
  
  if (vis->mode == 0 || vis->mode == 1)
    {
      rtRayTracing (rtAbsorptionCoefficients, net, vis);
    }
  else if (vis->mode == 2)
    {
      slStreamlines (slColourPalette, net, vis);
    }
  
  // here, intra-machine communications are handled through a binary
  // tree pattern and parallel pairwise blocking communications. The
  // master processor of the current machine gets the sub-images of
  // all the processors of that machine. Inter-machine communications,
  // needed if the number of machines is greater than one, take place
  // in the routine visRenderB.
  
  memcpy (vis->col_pixel_recv, vis->col_pixel_send, vis->col_pixels * sizeof(ColPixel));
  
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
	      net->err = MPI_Send (&vis->col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
	      if (vis->col_pixels > 0)
		{
#ifndef NOMPI
		  net->err = MPI_Send (&vis->col_pixel_send, vis->col_pixels, MPI_col_pixel_type,
					recv_id, 20, MPI_COMM_WORLD);
#endif
		}
	    }
	  else
	    {
#ifndef NOMPI
	      net->err = MPI_Recv (&col_pixels, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
				   net->status);
	      if (col_pixels > 0)
		{
		  net->err = MPI_Recv (&vis->col_pixel_send, col_pixels, MPI_col_pixel_type,
				       send_id, 20, MPI_COMM_WORLD, net->status);
		}
#else
	      col_pixels = 0;
#endif
	      for (n = 0; n < col_pixels; n++)
		{
		  col_pixel1 = &vis->col_pixel_send[ n ];
		  i = col_pixel1->i;
		  j = col_pixel1->j;
		  
		  if (*(col_pixel_id_p = &vis->col_pixel_id[ i * pixels_y + j ]) == -1)
		    {
		      if (vis->col_pixels == vis->col_pixels_max)
			{
			  vis->col_pixels_max <<= 1;
			  vis->col_pixel_recv = (ColPixel *)realloc(vis->col_pixel_recv,
								    sizeof(ColPixel) * vis->col_pixels_max);
			  vis->col_pixel_locked = (ColPixel *)realloc(vis->col_pixel_locked,
								      sizeof(ColPixel) * vis->col_pixels_max);
			}
		      col_pixel2 = &vis->col_pixel_recv[ *col_pixel_id_p = vis->col_pixels ];
		      
		      memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		      ++vis->col_pixels;
		    }
		  else
		    {
		      col_pixel2 = &vis->col_pixel_recv[ *col_pixel_id_p ];
		      
		      if (vis->mode == 0)
			{
			  col_pixel2->r += col_pixel1->r;
			  col_pixel2->g += col_pixel1->g;
			  col_pixel2->b += col_pixel1->b;
			}
		      else if (vis->mode == 1 || vis->mode == 2)
			{
			  if (col_pixel1->t < col_pixel2->t)
			    {
			      col_pixel2->r = col_pixel1->r;
			      col_pixel2->g = col_pixel1->g;
			      col_pixel2->b = col_pixel1->b;
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
		  memcpy (vis->col_pixel_send, vis->col_pixel_recv,
			  vis->col_pixels * sizeof(ColPixel));
		}
	    }
	  recv_id += comm_inc << 1;
	}
      comm_inc <<= 1;
    }
  
  if (net->machines == 1 || (net->id != 0 && net->id != master_proc_id))
    {
      return;
    }
  
  // inter-machine communications of sub-images begin here
  
  if (net->id != 0)
    {
      recv_id = 0;
#ifndef NOMPI
      net->err = MPI_Send (&vis->col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
      if (vis->col_pixels > 0)
	{
	  memcpy (vis->col_pixel_send, vis->col_pixel_recv,
		  vis->col_pixels * sizeof(ColPixel));
#ifndef NOMPI
	  net->err = MPI_Isend (vis->col_pixel_send,
				vis->col_pixels, MPI_col_pixel_type,
				recv_id, 30, MPI_COMM_WORLD,
				&net->req[ 1 ][ net->id * net->procs + recv_id ]);
#endif
	}
    }
  else
    {
      send_id = net->procs_per_machine[ net->id ];
      
      for (m = 1; m < net->machines; m++)
	{
#ifndef NOMPI
	  net->err = MPI_Recv (&vis->col_pixels_recv[ m-1 ], 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
			       net->status);
#endif
	  if (vis->col_pixels_recv[ m-1 ] > 0)
	    {
#ifndef NOMPI
	      net->err = MPI_Irecv (&vis->col_pixel_send[ (m-1) * (COLOURED_PIXELS_PER_PROC_MAX * sizeof(ColPixel)) ],
				    vis->col_pixels_recv[ m-1 ], MPI_col_pixel_type,
				    send_id, 30, MPI_COMM_WORLD,
				    &net->req[ 1 ][ (net->id + net->procs) * net->procs + send_id ]);
#endif
	    }
	  send_id += net->procs_per_machine[ m ];
	}
    }
}


void visRenderB (Net *net, Vis *vis)
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
  
  
  pixels_y = screen.pixels_y;
  
  if (net->machines > 1)
    {
      master_proc_id = 0;
      
      for (m = 0; m < net->machine_id[ net->id ]; m++)
	{
	  master_proc_id += net->procs_per_machine[ m ];
	}
      if (net->id != 0 && net->id != master_proc_id) return;
      
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
	  
	  for (m = 1; m < net->machines; m++)
	    {
	      if (vis->col_pixels_recv[ m-1 ] == 0)
		{
		  send_id += net->procs_per_machine[ m ];
		  continue;
		}
#ifndef NOMPI
	      net->err = MPI_Wait (&net->req[ 1 ][ (net->id + net->procs) * net->procs + send_id ], net->status);
#endif
	      offset = (m-1) * COLOURED_PIXELS_PER_PROC_MAX;
	      
	      for (n = 0; n < vis->col_pixels_recv[ m-1 ]; n++)
		{
		  col_pixel1 = &vis->col_pixel_send[ offset + n ];
		  i = col_pixel1->i;
		  j = col_pixel1->j;
		  
		  if (*(col_pixel_id_p = &vis->col_pixel_id[ i * pixels_y + j ]) == -1)
		    {
		      if (vis->col_pixels == vis->col_pixels_max)
			{
			  vis->col_pixels_max <<= 1;
			  vis->col_pixel_recv = (ColPixel *)realloc(vis->col_pixel_recv,
								    sizeof(ColPixel) * vis->col_pixels_max);
			}
		      col_pixel2 = &vis->col_pixel_recv[ *col_pixel_id_p = vis->col_pixels ];
		      
		      memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		      ++vis->col_pixels;
		    }
		  else
		    {
		      col_pixel2 = &vis->col_pixel_recv[ *col_pixel_id_p ];
		      
		      if (vis->mode == 0)
			{
			  col_pixel2->r += col_pixel1->r;
			  col_pixel2->g += col_pixel1->g;
			  col_pixel2->b += col_pixel1->b;
			}
		      else if (vis->mode == 1 || vis->mode == 2)
			{
			  if (col_pixel1->t < col_pixel2->t)
			    {
			      col_pixel2->r = col_pixel1->r;
			      col_pixel2->g = col_pixel1->g;
			      col_pixel2->b = col_pixel1->b;
			      col_pixel2->t = col_pixel1->t;
			    }
			}
		    }
		}
	    }
	}
    }
  
#ifndef BENCH
  
  float factor;
  
  if (net->id != 0) return;
  
  if (vis->mode == 0)
    {
      factor = 255.F * vis->absorption_factor;
    }
  else if (vis->mode == 1 || vis->mode == 2)
    {
      factor = 255.F;
    }
  
#ifdef RG
  
  memcpy (vis->col_pixel_locked, vis->col_pixel_recv,
  	  vis->col_pixels * sizeof(ColPixel));
  vis->col_pixels_locked = vis->col_pixels;
  
  for (n = 0; n < vis->col_pixels; n++)
    {
      vis->col_pixel_locked[ n ].r *= factor;
      vis->col_pixel_locked[ n ].g *= factor;
      vis->col_pixel_locked[ n ].b *= factor;
    }
  
#else
    
  FILE *image_file;
  XDR	xdr_image_file;
  
  short int pixel_i, pixel_j;
  
  unsigned char pixel_r, pixel_g, pixel_b;
  
  ColPixel *col_pixel_p;
  
  
  image_file = fopen (vis->image_file_name, "w");
  xdrstdio_create (&xdr_image_file, image_file, XDR_ENCODE);
  
  xdr_int (&xdr_image_file, &pixels_x);
  xdr_int (&xdr_image_file, &pixels_y);
  
  xdr_int (&xdr_image_file, &vis->col_pixels);
  
  for (n = 0; n < vis->col_pixels; n++)
    {
      col_pixel_p = &vis->col_pixel_recv[ n ];
      
      pixel_r = (unsigned char)max(0, min(255, (int)(255.F - factor * col_pixel_p->r)));
      pixel_g = (unsigned char)max(0, min(255, (int)(255.F - factor * col_pixel_p->g)));
      pixel_b = (unsigned char)max(0, min(255, (int)(255.F - factor * col_pixel_p->b)));
      
      pixel_i = col_pixel_p->i;
      pixel_j = col_pixel_p->j;
      
      xdr_u_char (&xdr_image_file, &pixel_r);
      xdr_u_char (&xdr_image_file, &pixel_g);
      xdr_u_char (&xdr_image_file, &pixel_b);
      xdr_short  (&xdr_image_file, &pixel_i);
      xdr_short  (&xdr_image_file, &pixel_j);
    }
  xdr_destroy (&xdr_image_file);
  fclose (image_file);
  
#endif // RG
#endif // BENCH
}


void visEnd (Net *net, Vis *vis)
{
  rtEnd (vis);
  
  slEnd (vis);
  
  /// free(vis->col_pixel_send);
  free(vis->col_pixel_id);
  free(vis->col_pixel_locked);
  free(vis->col_pixel_recv);
}

