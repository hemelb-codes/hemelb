#include "config.h"


void (*rtRayAABBIntersection[8]) (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t);

void (*rtTraverse[8]) (float org_x, float org_y, float org_z,
		       void (*rtAbsorptionCoefficients) (float flow_field_value, float t1, float t2,
							 float cutoff, float *r, float *g, float *b),
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


void rtTraverse000 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
	      	{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i -= vis->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	      
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= vis->block_size) < 0) return;
		  
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
    }
}

void rtTraverse001 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i      ) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j      ) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}	      
	      if ((i -= vis->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= vis->block_size) < 0) return;
	      
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse010 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i      ) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k      ) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i -= vis->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += vis->block_size) >= vis->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
    }
}

void rtTraverse011 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i      ) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i -= vis->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += vis->block_size) >= vis->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse100 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j      ) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k      ) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += vis->block_size2) >= vis->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= vis->block_size) < 0) return;
	      
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
    }
}

void rtTraverse101 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j      ) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += vis->block_size2) >= vis->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		      
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j -= vis->block_size) < 0) return;
	      
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse110 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k      ) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += vis->block_size2) >= vis->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		      
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += vis->block_size) >= vis->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (--k < 0) return;
	      
	      t = t_max_z;
	      t_max_z -= ray.inv_z;
	    }
	}
    }
}

void rtTraverse111 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, Vis *vis)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(vis->block_size_1, block_data.i));
  j = max(0, min(vis->block_size_1, block_data.j));
  k = max(0, min(vis->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= vis->shift2;
  j <<= vis->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + vis->flow_field_type ] *
	    vis->flow_field_value_max_inv[ vis->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_x, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((i += vis->block_size2) >= vis->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_y, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if ((j += vis->block_size) >= vis->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > vis->cutoff)
		{
		  AbsorptionCoefficients (value, vis->t_min = t, t_max_z, vis->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (vis->mode == 1) return;
		}
	      if (++k >= vis->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse000 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse000 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = block_data.min_x * ray.inv_x;
  t_max_y = block_data.min_y * ray.inv_y;
  t_max_z = block_data.min_z * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= vis->blocks_yz;
	      
	      block_data.min_x -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= vis->blocks_z;
	      
	      block_data.min_y -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse001 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse001 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x                  ) * ray.inv_x;
  t_max_y = (block_data.min_y                  ) * ray.inv_y;
  t_max_z = (block_data.min_z + vis->block_size) * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= vis->blocks_yz;
	      
	      block_data.min_x -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse001 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_x - block_data.min_z);
		  
		  rtTraverse001 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= vis->blocks_z;
	      
	      block_data.min_y -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse001 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse001 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
    }
}

void rtTraverse010 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse010 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x                  ) * ray.inv_x;
  t_max_y = (block_data.min_y + vis->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z                  ) * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= vis->blocks_yz;
	      
	      block_data.min_x -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse010 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse010 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += vis->blocks_z;
	      
	      block_data.min_y += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse010 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse010 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse011 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse011 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x                  ) * ray.inv_x;
  t_max_y = (block_data.min_y + vis->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z + vis->block_size) * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= vis->blocks_yz;
	      
	      block_data.min_x -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse011 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse011 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += vis->blocks_z;
	      
	      block_data.min_y += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse011 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse011 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
    }
}

void rtTraverse100 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse100 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + vis->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y                  ) * ray.inv_y;
  t_max_z = (block_data.min_z                  ) * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += vis->blocks_yz;
	      
	      block_data.min_x += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse100 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse100 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= vis->blocks_z;
	      
	      block_data.min_y -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse100 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse100 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse101 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse101 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + vis->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y                  ) * ray.inv_y;
  t_max_z = (block_data.min_z + vis->block_size) * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += vis->blocks_yz;
	      
	      block_data.min_x += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse101 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse101 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= vis->blocks_z;
	      
	      block_data.min_y -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse101 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse101 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
    }
}

void rtTraverse110 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse110 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + vis->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y + vis->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z                  ) * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += vis->blocks_yz;
	      
	      block_data.min_x += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse110 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse110 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += vis->blocks_z;
	      
	      block_data.min_y += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse110 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse110 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse111 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, Vis *vis)
{
  float dx, dy, dz;
  float t_max_x, t_max_y, t_max_z;
  float t_delta_x, t_delta_y, t_delta_z;
  
  int i, j, k;
  int ii, jj, kk;
  
  BlockData block_data;
  
  
  dx = org_x - cluster_p->block_min_x;
  dy = org_y - cluster_p->block_min_y;
  dz = org_z - cluster_p->block_min_z;
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(vis->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(vis->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(vis->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * vis->block_size) - dx;
  block_data.min_y = (float)(j * vis->block_size) - dy;
  block_data.min_z = (float)(k * vis->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= vis->blocks_yz;
  jj *= vis->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse111 (block_data, AbsorptionCoefficients, net, vis);
      
      if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + vis->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y + vis->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z + vis->block_size) * ray.inv_z;
  
  t_delta_x = vis->block_size * ray.inv_x;
  t_delta_y = vis->block_size * ray.inv_y;
  t_delta_z = vis->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += vis->blocks_yz;
	      
	      block_data.min_x += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse111 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse111 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += vis->blocks_z;
	      
	      block_data.min_y += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse111 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += vis->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse111 (block_data, AbsorptionCoefficients, net, vis);
		  
		  if (vis->mode == 1 && vis->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
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


void rtRayTracing (void (*AbsorptionCoefficients) (float flow_field_data, float t1, float t2,
						   float cutoff, float *r, float *g, float *b),
		   Net *net, Vis *vis)
{
  // here, the ray tracing is performed and the intra-machine communications take place
  
  float screen_vtx_x, screen_vtx_y, screen_vtx_z;
  
  float px0, py0, pz0;
  float px1, py1, pz1;
  float nx0, ny0, nz0;
  float t;
  float temp1;
  
  float par1x, par1y, par1z;
  float par2x, par2y, par2z;
  float par3x, par3y, par3z;
  
  float v_min_x, v_min_y;
  float v_max_x, v_max_y;
  float vx[2], vy[2], vz[2];
  float scale_x, scale_y;
  
  int pixels_x, pixels_y;
  int i, j, k;
  int cluster_id;
  int min_i, min_j, max_i, max_j;
  int viewpoint_flag;
  int ray_dir_code;
  
  AABB aabb;
  
  Cluster *cluster_p;
  
  
  px0 = viewpoint.pos_x;
  py0 = viewpoint.pos_y;
  pz0 = viewpoint.pos_z;
  
  screen_vtx_x = EPSILON + screen.ctr_x - screen.dir1x - screen.dir2x - px0;
  screen_vtx_y = EPSILON + screen.ctr_y - screen.dir1y - screen.dir2y - py0;
  screen_vtx_z = EPSILON + screen.ctr_z - screen.dir1z - screen.dir2z - pz0;
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  par1x = screen.dir1x * (2.F / (float)pixels_x);
  par1y = screen.dir1y * (2.F / (float)pixels_x);
  par1z = screen.dir1z * (2.F / (float)pixels_x);
  
  par2x = screen.dir2x * (2.F / (float)pixels_y);
  par2y = screen.dir2y * (2.F / (float)pixels_y);
  par2z = screen.dir2z * (2.F / (float)pixels_y);
  
  scale_x = 1.F / screen.par_x;
  scale_y = 1.F / screen.par_y;

  for (cluster_id = 0; cluster_id < vis->clusters; cluster_id++)
    {
      cluster_p = &vis->cluster[ cluster_id ];
      
      // the image-based projection of the cluster bounding box is calculated here
      
      vx[0] = cluster_p->block_min_x;
      vy[0] = cluster_p->block_min_y;
      vz[0] = cluster_p->block_min_z;
      vx[1] = vx[0] + cluster_p->blocks_x * vis->block_size;
      vy[1] = vy[0] + cluster_p->blocks_y * vis->block_size;
      vz[1] = vz[0] + cluster_p->blocks_z * vis->block_size;
      
      v_min_x = 1.e+30F;
      v_min_y = 1.e+30F;
      v_max_x = -1.e+30F;
      v_max_y = -1.e+30F;
      
      for (i = 0; i < 2; i++)
	{
	  for (j = 0; j < 2; j++)
	    {
	      for (k = 0; k < 2; k++)
		{
		  visProject (vx[i], vy[j], vz[k], &px1, &py1, &pz1);
		  
		  v_min_x = fminf(v_min_x, px1);
		  v_min_y = fminf(v_min_y, py1);
		  v_max_x = fmaxf(v_max_x, px1);
		  v_max_y = fmaxf(v_max_y, py1);
		}
	    }
	}
      min_i = (int)(scale_x * (v_min_x + screen.max_x));
      max_i = (int)(scale_x * (v_max_x + screen.max_x));
      min_j = (int)(scale_y * (v_min_y + screen.max_y));
      max_j = (int)(scale_y * (v_max_y + screen.max_y));
      
      if (min_i >= pixels_x || min_j >= pixels_y || max_i < 0 || max_j < 0) continue;
      
      min_i = max(min_i, 0);
      min_j = max(min_j, 0);
      max_i = min(max_i, pixels_x - 1);
      max_j = min(max_j, pixels_y - 1);
      
      if (px0 >= vx[0] && py0 >= vy[0] && pz0 >= vz[0] &&
	  px0 <= vx[1] && py0 <= vy[1] && pz0 <= vz[1])
	{
	  viewpoint_flag = 1;
	}
      else
	{
	  viewpoint_flag = 0;
	}
      aabb.acc_1 = vx[1] - px0;
      aabb.acc_2 = vx[0] - px0;
      aabb.acc_3 = vy[1] - py0;
      aabb.acc_4 = vy[0] - py0;
      aabb.acc_5 = vz[1] - pz0;
      aabb.acc_6 = vz[0] - pz0;
      
      par3x = screen_vtx_x + (float)(min_i + 0.5F) * par1x + (float)(min_j + 0.5F) * par2x;
      par3y = screen_vtx_y + (float)(min_i + 0.5F) * par1y + (float)(min_j + 0.5F) * par2y;
      par3z = screen_vtx_z + (float)(min_i + 0.5F) * par1z + (float)(min_j + 0.5F) * par2z;
      
      for (i = min_i; i <= max_i; i++)
	{
	  nx0 = par3x;
	  ny0 = par3y;
	  nz0 = par3z;
	  
	  for (j = min_j; j <= max_j; j++)
	    {
	      px1 = px0;
	      py1 = py0;
	      pz1 = pz0;
	      
	      ray.dir_x = nx0;
	      ray.dir_y = ny0;
	      ray.dir_z = nz0;
	      temp1 = 1.F / sqrtf(nx0 * nx0 + ny0 * ny0 + nz0 * nz0);
	      ray.dir_x *= temp1;
	      ray.dir_y *= temp1;
	      ray.dir_z *= temp1;
	      
	      ray.inv_x = 1.F / ray.dir_x;
	      ray.inv_y = 1.F / ray.dir_y;
	      ray.inv_z = 1.F / ray.dir_z;
	      
	      ray_dir_code = ray.dir_z > 0.F;
	      ray_dir_code |= (ray.dir_y > 0.F) << 1;
	      ray_dir_code |= (ray.dir_x > 0.F) << 2;
	      
	      nx0 += par2x;
	      ny0 += par2y;
	      nz0 += par2z;
	      
	      t = 0.F;
	      
	      if (!viewpoint_flag)
		{
		  (*rtRayAABBIntersection[ ray_dir_code ]) (&aabb, ray.inv_x, ray.inv_y, ray.inv_z, &t);
		  
		  if (t >= 1.e+30F) continue;
		  
		  px1 += t * ray.dir_x;
		  py1 += t * ray.dir_y;
		  pz1 += t * ray.dir_z;
		}
	      
	      ray.col_r = 0.F;
	      ray.col_g = 0.F;
	      ray.col_b = 0.F;
	      
	      vis->t_min = 1.e+30F;
	      
	      (*rtTraverse[ ray_dir_code ]) (px1, py1, pz1, AbsorptionCoefficients, cluster_p, net, vis);
	      
	      if (vis->t_min >= 1.e+30F) continue;
	      
	      visWritePixel (ray.col_r, ray.col_g, ray.col_b, vis->t_min + t, i, j, vis);
	    }
	  par3x += par1x;
	  par3y += par1y;
	  par3z += par1z;
	}
    }
}


void rtEnd (Vis *vis)
{
  free(vis->cluster);
}


void slInit (Net *net, Vis *vis)
{
  float *seed_p;
  
  int boundary_sites_count;
  int site_i, site_j, site_k;
  int i, j, k;
  int m, n;
  
  unsigned int site_id, site_type;
  
  DataBlock *map_block_p;
  
  
  vis->seeds_max = 10000;
  vis->seed = (float *)malloc(sizeof(float) * 3 * vis->seeds_max);
  
  vis->seeds = 0;
  
  boundary_sites_count = 0;
  
  n = -1;
  
  for (i = 0; i < net->sites_x; i += net->block_size)
    {
      for (j = 0; j < net->sites_y; j += net->block_size)
	{
	  for (k = 0; k < net->sites_z; k += net->block_size)
	    {
	      if (net->proc_id[ ++n ] != net->id)
		{
		  continue;
		}
	      map_block_p = &net->map_block[ n ];
	      
	      m = -1;
	      
	      for (site_i = i; site_i < i + net->block_size; site_i++)
		{
		  for (site_j = j; site_j < j + net->block_size; site_j++)
		    {
		      for (site_k = k; site_k < k + net->block_size; site_k++)
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
			  seed_p = &vis->seed[ 3 * vis->seeds ];
			  seed_p[ 0 ] = (float)site_i + (0.5F - vis->half_dim_x);
			  seed_p[ 1 ] = (float)site_j + (0.5F - vis->half_dim_y);
			  seed_p[ 2 ] = (float)site_k + (0.5F - vis->half_dim_z);
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


void slStreamlines (void ColourPalette (float vel_m, float *r, float *g, float *b),
		    Net *net, Vis *vis)
{
  double density, vx, vy, vz;
  
  float scale_x, scale_y;
  float vel_m;
  float r, g, b;
  float z_old, z_new;
  float x, y, z;
  float z1;
  float z2;
  float dz;
  float d, incE, incNE;
  float *streamline_p;
  
  int pixels_x, pixels_y;
  int cycles;
  int is_inside_my_subdomain;
  int proc_id, neigh_proc_index;
  int site_i, site_j, site_k;
  int bi, bj, bk;
  int si, sj, sk;
  int block_id;
  int i_old, j_old;
  int i_new, j_new;
  int i1, j1, i2, j2;
  int di, dj;
  int i, j;
  int m, n;
  int streamlines_max;
  
  unsigned int site_id;
  
  DataBlock *map_block_p;
  
  NeighProc *neigh_proc_p;
  
  
  vis->streamlines = vis->seeds;
  memcpy (vis->streamline, vis->seed, vis->seeds * 3 * sizeof(float));
  
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
	  streamline_p = &vis->streamline[ 3 * n ];
	  
	  site_i = (int)(streamline_p[ 0 ] + vis->half_dim_x);
	  site_j = (int)(streamline_p[ 1 ] + vis->half_dim_y);
	  site_k = (int)(streamline_p[ 2 ] + vis->half_dim_z);
	  
	  if (site_i < 0 || site_i >= net->sites_x ||
	      site_j < 0 || site_j >= net->sites_y ||
	      site_k < 0 || site_k >= net->sites_z)
	    {
	      --vis->streamlines;
	      continue;
	    }
	  bi = site_i >> net->shift;
	  bj = site_j >> net->shift;
	  bk = site_k >> net->shift;
	  
	  block_id = (bi * net->blocks_y + bj) * net->blocks_z + bk;
	  
	  proc_id = net->proc_id[ block_id ];
	  
	  if (proc_id != net->id)
	    {
	      --vis->streamlines;
	      continue;
	    }
	  si = site_i - (bi << net->shift);
	  sj = site_j - (bj << net->shift);
	  sk = site_k - (bk << net->shift);
	  
	  site_id = net->map_block[ block_id ].site_data[ (((si << net->shift) + sj) << net->shift) + sk ];
	  
	  if (site_id & (1U << 31U))
	    {
	      --vis->streamlines;
	      continue;
	    }
	  visProject (streamline_p[ 0 ], streamline_p[ 1 ], streamline_p[ 2 ], &x, &y, &z);
	  
	  i_old = (int)(scale_x * (x + screen.max_x));
	  j_old = (int)(scale_y * (y + screen.max_y));
	  z_old = z;
	  
	  is_inside_my_subdomain = 1;
	  
	  while (is_inside_my_subdomain)
	    {
	      is_inside_my_subdomain = 0;
	      
	      lbmDensityAndVelocity (&f_old[ site_id*15 ], &density, &vx, &vy, &vz);
	      
	      vel_m = (float)sqrt(vx * vx + vy * vy + vz * vz);
	      
	      if (vel_m < 1.e-6)
		{
		  --vis->streamlines;
		  continue;
		}
	      ColourPalette (vel_m * vis->flow_field_value_max_inv[ VELOCITY ],
			     &r, &g, &b);
	      
	      vel_m = 1.F / vel_m;
	      streamline_p[ 0 ] += (float)vx * vel_m;
	      streamline_p[ 1 ] += (float)vy * vel_m;
	      streamline_p[ 2 ] += (float)vz * vel_m;
	      
	      visProject (streamline_p[ 0 ], streamline_p[ 1 ], streamline_p[ 2 ], &x, &y, &z);
	      
	      i_new = (int)(scale_x * (x + screen.max_x));
	      j_new = (int)(scale_y * (y + screen.max_y));
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
			  visWritePixel (r, g, b, z, i, j, vis);
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
			  visWritePixel (r, g, b, z, i, j, vis);
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
	      site_i = (int)(streamline_p[ 0 ] + vis->half_dim_x);
	      site_j = (int)(streamline_p[ 1 ] + vis->half_dim_y);
	      site_k = (int)(streamline_p[ 2 ] + vis->half_dim_z);
	      
	      if (site_i < 0 || site_i >= net->sites_x ||
		  site_j < 0 || site_j >= net->sites_y ||
		  site_k < 0 || site_k >= net->sites_z)
		{
		  --vis->streamlines;
		  continue;
		}
	      bi = site_i >> net->shift;
	      bj = site_j >> net->shift;
	      bk = site_k >> net->shift;
	      
	      block_id = (bi * net->blocks_y + bj) * net->blocks_z + bk;
	      
	      proc_id = net->proc_id[ block_id ];
	      
	      if (proc_id != net->id && proc_id != (1 << 14))
		{
		  neigh_proc_index = net->from_proc_id_to_neigh_proc_index[ proc_id ];
		  
		  if (streamlines_to_send[ neigh_proc_index ] < STREAMLINES_MAX)
		    {
		      memcpy (&streamline_to_send[ neigh_proc_index ][ 3*streamlines_to_send[neigh_proc_index] ],
		  	      streamline_p, 3 * sizeof(float));
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
	      si = site_i - (bi << net->shift);
	      sj = site_j - (bj << net->shift);
	      sk = site_k - (bk << net->shift);
	      
	      site_id = map_block_p->site_data[ (((si << net->shift) + sj) << net->shift) + sk ];
	      
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
      for (m = 0; m < net->neigh_procs; m++)
      	{
      	  neigh_proc_p = &net->neigh_proc[ m ];
      	  
      	  net->err = MPI_Ssend (&streamlines_to_send[m], 1, MPI_INT,
				neigh_proc_p->id, 30, MPI_COMM_WORLD);
      	  net->err = MPI_Recv (&streamlines_to_recv[m], 1, MPI_INT,
			       neigh_proc_p->id, 30, MPI_COMM_WORLD, net->status);
      	  
	  if (streamlines_to_send[ m ] > 0)
      	    {
      	      net->err = MPI_Ssend (&streamline_to_send[m][0], 3 * streamlines_to_send[ m ],
				    MPI_REAL, neigh_proc_p->id, 30, MPI_COMM_WORLD);
      	    }
	  if (streamlines_to_recv[ m ] > 0)
      	    {
      	      net->err = MPI_Recv (&streamline_to_recv[m][0], 3 * streamlines_to_recv[ m ],
				   MPI_REAL, neigh_proc_p->id, 30, MPI_COMM_WORLD, net->status);
	      
              streamlines_max = min(STREAMLINES_MAX - vis->streamlines,
      	  			    streamlines_to_recv[ m ]);
      	      
      	      memcpy (&vis->streamline[ vis->streamlines ], &streamline_to_recv[ m ][ 0 ],
      	  	      streamlines_max * 3 * sizeof(float));
      	      vis->streamlines += streamlines_max;
      	    }
	}
      
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


void visProject (float px1, float py1, float pz1, float *px2, float *py2, float *pz2)
{
  float x1, y1, z1, x2, y2, z2;
  float temp;
  
  
  x1 = px1 - viewpoint.pos_x;
  y1 = py1 - viewpoint.pos_y;
  z1 = pz1 - viewpoint.pos_z;
  
  temp = viewpoint.cos_1 * z1 + viewpoint.sin_1 * x1;
  
  x2 = viewpoint.cos_1 * x1 - viewpoint.sin_1 * z1;
  y2 = viewpoint.cos_2 * y1 - viewpoint.sin_2 * temp;
  z2 = viewpoint.cos_2 * temp + viewpoint.sin_2 * y1;
  
  *pz2 = -z2;
  
  temp = screen.dist / *pz2;
  
  *px2 = temp * x2;
  *py2 = temp * y2;
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
  
  viewpoint.pos_x = temp * viewpoint.sin_1 + ctr_x;
  viewpoint.pos_y = rad  * viewpoint.sin_2 + ctr_y;
  viewpoint.pos_z = temp * viewpoint.cos_1 + ctr_z;
  
  temp = dist / rad;
  
  screen.ctr_x = viewpoint.pos_x + temp * (ctr_x - viewpoint.pos_x);
  screen.ctr_y = viewpoint.pos_y + temp * (ctr_y - viewpoint.pos_y);
  screen.ctr_z = viewpoint.pos_z + temp * (ctr_z - viewpoint.pos_z);
  
  screen.zoom = zoom;
  screen.dist = dist;
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     screen.max_x, 0.0F, 0.0F,
	     &screen.dir1x, &screen.dir1y, &screen.dir1z);
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     0.0F, screen.max_y, 0.0F,
	     &screen.dir2x, &screen.dir2y, &screen.dir2z);
  
  screen.par_x = (2.F * screen.max_x) / (float)pixels_x;
  screen.par_y = (2.F * screen.max_y) / (float)pixels_y;
}


void visWritePixel (float r, float g, float b, float t, int i, int j, Vis *vis)
{
  int pixels_x, pixels_y;
  int *col_pixel_id_p;
  int err;
  
  ColPixel *col_pixel_p;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  if (*(col_pixel_id_p = &vis->col_pixel_id[ i * pixels_y + j ]) == -1)
    {
      if (vis->col_pixels >= COLOURED_PIXELS_PER_PROC_MAX)
	{
	  printf (" too many coloured pixels per proc\n");
	  printf (" the execution is terminated\n");
	  err = MPI_Abort (MPI_COMM_WORLD, 1);
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
      col_pixel_p->r += r;
      col_pixel_p->g += g;
      col_pixel_p->b += b;
    }
  else if (vis->mode == 1 || vis->mode == 2)
    {
      if (t < col_pixel_p->t)
	{
	  col_pixel_p->r = r;
	  col_pixel_p->g = g;
	  col_pixel_p->b = b;
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
      parameters_file = fopen (parameters_file_name, "r");
      
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
  net->err = MPI_Bcast (par_to_send, 14, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
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
  int col_pixel_count = 6;
  int col_pixel_blocklengths[6] = {1, 1, 1, 1, 1, 1};
  
  MPI_Datatype col_pixel_types[6] = {MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL,
				     MPI_SHORT, MPI_SHORT};
  
  MPI_Aint col_pixel_disps[6];
  
  
  vis->image_file_name = image_file_name;
  
  vis->col_pixels_max = 512 * 512;
  
  /// vis->col_pixel_send = (ColPixel *)malloc(sizeof(ColPixel) * vis->col_pixels_max);
  
  vis->col_pixel_recv = (ColPixel *)malloc(sizeof(ColPixel) * vis->col_pixels_max);
  vis->col_pixel_locked = (ColPixel *)malloc(sizeof(ColPixel) * vis->col_pixels_max);
  
  vis->pixels_max = IMAGE_SIZE;
  vis->col_pixel_id = (int *)malloc(sizeof(int) * vis->pixels_max);
  
  // create the derived datatype for the MPI communications
  
  col_pixel_disps[0] = 0;
  
  for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_INTEGER)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
	}
      else if (col_pixel_types[i - 1] == MPI_DOUBLE)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(double) * col_pixel_blocklengths[i - 1]);
	}
      else if (col_pixel_types[i - 1] == MPI_REAL)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(float) * col_pixel_blocklengths[i - 1]);
	}
      else if (col_pixel_types[i - 1] == MPI_SHORT)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(short int) * col_pixel_blocklengths[i - 1]);
	}
    }
  MPI_Type_struct (col_pixel_count, col_pixel_blocklengths, col_pixel_disps, col_pixel_types, &MPI_col_pixel_type);
  MPI_Type_commit (&MPI_col_pixel_type);
  
  
  rtInit (net, vis);
  
  slInit (net, vis);
}


void visRenderA (void (*rtAbsorptionCoefficients) (float flow_field_data, float t1, float t2,
						   float cutoff, float *r, float *g, float *b),
		 void slColourPalette (float vel_m, float *r, float *g, float *b),
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
	      net->err = MPI_Ssend (&vis->col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
	      
	      if (vis->col_pixels > 0)
		{
		  net->err = MPI_Ssend (&vis->col_pixel_send, vis->col_pixels, MPI_col_pixel_type, recv_id, 20, MPI_COMM_WORLD);
		}
	    }
	  else
	    {
	      net->err = MPI_Recv (&col_pixels, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
				   net->status);
	      
	      if (col_pixels > 0)
		{
		  net->err = MPI_Recv (&vis->col_pixel_send, col_pixels, MPI_col_pixel_type, send_id, 20, MPI_COMM_WORLD,
				       net->status);
		}
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
      
      net->err = MPI_Ssend (&vis->col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
      
      if (vis->col_pixels > 0)
	{
	  memcpy (vis->col_pixel_send, vis->col_pixel_recv,
		  vis->col_pixels * sizeof(ColPixel));
	  
	  net->err = MPI_Isend (vis->col_pixel_send,
				vis->col_pixels, MPI_col_pixel_type,
				recv_id, 30, MPI_COMM_WORLD,
				&net->req[ 1 ][ net->id * net->procs + recv_id ]);
	}
    }
  else
    {
      send_id = net->procs_per_machine[ net->id ];
      
      for (m = 1; m < net->machines; m++)
	{
	  net->err = MPI_Recv (&vis->col_pixels_recv[ m-1 ], 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
			       net->status);
	  
	  if (vis->col_pixels_recv[ m-1 ] > 0)
	    {
	      net->err = MPI_Irecv (&vis->col_pixel_send[ (m-1) * (COLOURED_PIXELS_PER_PROC_MAX * sizeof(ColPixel)) ],
				    vis->col_pixels_recv[ m-1 ], MPI_col_pixel_type,
				    send_id, 30, MPI_COMM_WORLD,
				    &net->req[ 1 ][ (net->id + net->procs) * net->procs + send_id ]);
	    }
	  send_id += net->procs_per_machine[ m ];
	}
    }
}


void visRenderB (Net *net, Vis *vis)
{
  // here, the intra-machine communications take place and the buffer
  // to stream to the client or the output image are set
  
  float factor;
  
  int pixels_x, pixels_y;
  int i, j;
  int m, n;
  int *col_pixel_id_p;
  int send_id, recv_id;
  int master_proc_id;
  int offset;

  ColPixel *col_pixel1, *col_pixel2;
  
  
  pixels_x = screen.pixels_x;
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
	  
	  net->err = MPI_Wait (&net->req[ 1 ][ net->id * net->procs + recv_id ], net->status);
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
	      net->err = MPI_Wait (&net->req[ 1 ][ (net->id + net->procs) * net->procs + send_id ], net->status);
	      
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
