#include "config.h"


void (*rtRayAABBIntersection[8]) (AABB *aabb, float inv_x, float inv_y, float inv_z, float *t);

void (*rtTraverse[8]) (float org_x, float org_y, float org_z,
		       void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						       float cutoff, float *r, float *g, float *b),
		       Cluster *cluster_p, Net *net, RT *rt);


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
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
	      	{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((i -= rt->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	      
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
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
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j -= rt->block_size) < 0) return;
		  
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
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
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i      ) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j      ) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}	      
	      if ((i -= rt->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j -= rt->block_size) < 0) return;
	      
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse010 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i      ) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k      ) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((i -= rt->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
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
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j += rt->block_size) >= rt->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
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
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)i      ) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((i -= rt->block_size2) < 0) return;
	      
	      t = t_max_x;
	      t_max_x -= ray.inv_x;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j += rt->block_size) >= rt->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse100 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j      ) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k      ) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((i += rt->block_size2) >= rt->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
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
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j -= rt->block_size) < 0) return;
	      
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
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
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)j      ) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((i += rt->block_size2) >= rt->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		      
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j -= rt->block_size) < 0) return;
	      
	      t = t_max_y;
	      t_max_y -= ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse110 (BlockData block_data,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)k      ) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((i += rt->block_size2) >= rt->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		      
		  if (rt->is_isosurface) return;
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
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j += rt->block_size) >= rt->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
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
						    float cutoff, float *r, float *g, float *b), Net *net, RT *rt)
{
  float t_max_x, t_max_y, t_max_z;
  float t;
  float value;
  
  int i, j, k;
  
  unsigned int site_data;
  
  
  t = block_data.t;
  
  i = max(0, min(rt->block_size_1, block_data.i));
  j = max(0, min(rt->block_size_1, block_data.j));
  k = max(0, min(rt->block_size_1, block_data.k));
  
  t_max_x = (block_data.min_x + (float)(i + 1)) * ray.inv_x;
  t_max_y = (block_data.min_y + (float)(j + 1)) * ray.inv_y;
  t_max_z = (block_data.min_z + (float)(k + 1)) * ray.inv_z;
  
  i <<= rt->shift2;
  j <<= rt->shift;
  
  for (;;)
    {
      value = -1.F;
      
      if ((site_data = block_data.p->site_data[ i + j + k ]) != (1U << 31U))
	{
	  value = flow_field[ 3 * site_data + rt->flow_field_type ] *
	    rt->flow_field_value_max_inv[ rt->flow_field_type ];
	}
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_x, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((i += rt->block_size2) >= rt->block_size3) return;
	      
	      t = t_max_x;
	      t_max_x += ray.inv_x;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_y, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if ((j += rt->block_size) >= rt->block_size2) return;
	      
	      t = t_max_y;
	      t_max_y += ray.inv_y;
	    }
	  else
	    {
	      if (value > rt->cutoff)
		{
		  AbsorptionCoefficients (value, rt->t_min = t, t_max_z, rt->cutoff, &ray.col_r, &ray.col_g, &ray.col_b);
		  
		  if (rt->is_isosurface) return;
		}
	      if (++k >= rt->block_size) return;
	      
	      t = t_max_z;
	      t_max_z += ray.inv_z;
	    }
	}
    }
}

void rtTraverse000 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse000 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = block_data.min_x * ray.inv_x;
  t_max_y = block_data.min_y * ray.inv_y;
  t_max_z = block_data.min_z * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= rt->blocks_yz;
	      
	      block_data.min_x -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= rt->blocks_z;
	      
	      block_data.min_y -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse000 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse001 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse001 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x                 ) * ray.inv_x;
  t_max_y = (block_data.min_y                 ) * ray.inv_y;
  t_max_z = (block_data.min_z + rt->block_size) * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= rt->blocks_yz;
	      
	      block_data.min_x -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse001 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_x - block_data.min_z);
		  
		  rtTraverse001 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= rt->blocks_z;
	      
	      block_data.min_y -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse001 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse001 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
    }
}

void rtTraverse010 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse010 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x                 ) * ray.inv_x;
  t_max_y = (block_data.min_y + rt->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z                 ) * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= rt->blocks_yz;
	      
	      block_data.min_x -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse010 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse010 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += rt->blocks_z;
	      
	      block_data.min_y += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse010 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse010 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse011 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse011 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x                 ) * ray.inv_x;
  t_max_y = (block_data.min_y + rt->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z + rt->block_size) * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i -= cluster_p->blocks_yz) < 0) return;
	      ii -= rt->blocks_yz;
	      
	      block_data.min_x -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse011 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x -= t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse011 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += rt->blocks_z;
	      
	      block_data.min_y += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse011 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse011 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
    }
}

void rtTraverse100 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse100 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + rt->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y                 ) * ray.inv_y;
  t_max_z = (block_data.min_z                 ) * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += rt->blocks_yz;
	      
	      block_data.min_x += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse100 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse100 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= rt->blocks_z;
	      
	      block_data.min_y -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse100 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse100 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse101 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse101 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + rt->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y                 ) * ray.inv_y;
  t_max_z = (block_data.min_z + rt->block_size) * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += rt->blocks_yz;
	      
	      block_data.min_x += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse101 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse101 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j -= cluster_p->blocks_z) < 0) return;
	      jj -= rt->blocks_z;
	      
	      block_data.min_y -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse101 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y -= t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse101 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
    }
}

void rtTraverse110 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse110 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + rt->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y + rt->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z                 ) * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += rt->blocks_yz;
	      
	      block_data.min_x += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse110 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse110 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += rt->blocks_z;
	      
	      block_data.min_y += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse110 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (--k < 0) return;
	      --kk;
	      
	      block_data.min_z -= rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse110 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z -= t_delta_z;
	    }
	}
    }
}

void rtTraverse111 (float org_x, float org_y, float org_z,
		    void (*AbsorptionCoefficients) (float flow_field_value, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		    Cluster *cluster_p, Net *net, RT *rt)
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
  
  i = max(0, min(cluster_p->blocks_x - 1, (int)(rt->block_size_inv * dx)));
  j = max(0, min(cluster_p->blocks_y - 1, (int)(rt->block_size_inv * dy)));
  k = max(0, min(cluster_p->blocks_z - 1, (int)(rt->block_size_inv * dz)));
  
  block_data.min_x = (float)(i * rt->block_size) - dx;
  block_data.min_y = (float)(j * rt->block_size) - dy;
  block_data.min_z = (float)(k * rt->block_size) - dz;
  
  ii = i + cluster_p->block_min_i;
  jj = j + cluster_p->block_min_j;
  kk = k + cluster_p->block_min_k;
  
  ii *= rt->blocks_yz;
  jj *= rt->blocks_z;
  
  if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
    {
      block_data.t = 0.F;
      block_data.i = (int)(-block_data.min_x);
      block_data.j = (int)(-block_data.min_y);
      block_data.k = (int)(-block_data.min_z);
      
      rtTraverse111 (block_data, AbsorptionCoefficients, net, rt);
      
      if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
    }
  i *= cluster_p->blocks_yz;
  j *= cluster_p->blocks_z;
  
  t_max_x = (block_data.min_x + rt->block_size) * ray.inv_x;
  t_max_y = (block_data.min_y + rt->block_size) * ray.inv_y;
  t_max_z = (block_data.min_z + rt->block_size) * ray.inv_z;
  
  t_delta_x = rt->block_size * ray.inv_x;
  t_delta_y = rt->block_size * ray.inv_y;
  t_delta_z = rt->block_size * ray.inv_z;
  
  for (;;)
    {
      if (t_max_x < t_max_y)
	{
	  if (t_max_x < t_max_z)
	    {
	      if ((i += cluster_p->blocks_yz) >= cluster_p->blocks) return;
	      ii += rt->blocks_yz;
	      
	      block_data.min_x += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_x;
		  block_data.i = (int)(t_max_x * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_x * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_x * ray.dir_z - block_data.min_z);

		  rtTraverse111 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_x += t_delta_x;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse111 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
      else
	{
	  if (t_max_y < t_max_z)
	    {
	      if ((j += cluster_p->blocks_z) >= cluster_p->blocks_yz) return;
	      jj += rt->blocks_z;
	      
	      block_data.min_y += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_y;
		  block_data.i = (int)(t_max_y * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_y * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_y * ray.dir_z - block_data.min_z);
		  
		  rtTraverse111 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_y += t_delta_y;
	    }
	  else
	    {
	      if (++k >= cluster_p->blocks_z) return;
	      ++kk;
	      
	      block_data.min_z += rt->block_size;
	      
	      if ((block_data.p = &net->map_block[ ii + jj + kk ])->site_data != NULL)
		{
		  block_data.t = t_max_z;
		  block_data.i = (int)(t_max_z * ray.dir_x - block_data.min_x);
		  block_data.j = (int)(t_max_z * ray.dir_y - block_data.min_y);
		  block_data.k = (int)(t_max_z * ray.dir_z - block_data.min_z);
		  
		  rtTraverse111 (block_data, AbsorptionCoefficients, net, rt);
		  
		  if (rt->is_isosurface && rt->t_min < 1.e+30F) return;
		}
	      t_max_z += t_delta_z;
	    }
	}
    }
}


void rtProject (float px1, float py1, float pz1, float *px2, float *py2)
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
  
  temp = screen.dist / -z2;
  
  *px2 = temp * x2;
  *py2 = temp * y2;
}


void rtRayTracingA (void (*AbsorptionCoefficients) (float flow_field_data, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		   Net *net, RT *rt)
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
  int m, n;
  int cluster_id;
  int min_i, min_j, max_i, max_j;
  int viewpoint_flag;
  int ray_dir_code;
  int *col_pixel_id_p;
  int col_pixels;
  int comm_inc, send_id, recv_id;
  int machine_id, master_proc_id;
  
  AABB aabb;
  
  ColPixel *col_pixel_p, *col_pixel1, *col_pixel2;
  
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
  
  if (pixels_x * pixels_y > rt->pixels_max)
    {
      rt->pixels_max = pixels_x * pixels_y;
      
      rt->col_pixel_id = (int *)realloc(rt->col_pixel_id, sizeof(int) * rt->pixels_max);
      
#ifndef BENCH
#ifdef RG
      if (net->id == 0)
	{
	  pixel_data = (unsigned char *)realloc(pixel_data, sizeof(unsigned char) * 3 * pixels_x * pixels_y);
	  compressed_data = (unsigned char *)realloc(compressed_data, sizeof(unsigned char) * 3 * pixels_x * pixels_y);
	}
#endif // RG
#endif // BENCH
    }

  for (i = 0; i < pixels_x * pixels_y; i++)
    {
      rt->col_pixel_id[ i ] = -1;
    }
  rt->col_pixels = 0;

  for (cluster_id = 0; cluster_id < rt->clusters; cluster_id++)
    {
      cluster_p = &rt->cluster[ cluster_id ];
      
      // the image-based projection of the cluster bounding box is calculated here
      
      vx[0] = cluster_p->block_min_x;
      vy[0] = cluster_p->block_min_y;
      vz[0] = cluster_p->block_min_z;
      vx[1] = vx[0] + cluster_p->blocks_x * rt->block_size;
      vy[1] = vy[0] + cluster_p->blocks_y * rt->block_size;
      vz[1] = vz[0] + cluster_p->blocks_z * rt->block_size;
      
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
		  rtProject (vx[i], vy[j], vz[k], &px1, &py1);
		  
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
	      
	      rt->t_min = 1.e+30F;
	      
	      (*rtTraverse[ ray_dir_code ]) (px1, py1, pz1, AbsorptionCoefficients, cluster_p, net, rt);

	      if (rt->t_min >= 1.e+30F) continue;

	      if (*(col_pixel_id_p = &rt->col_pixel_id[ i * pixels_y + j ]) == -1)
		{
		  if (rt->col_pixels >= COLOURED_PIXELS_PER_PROC_MAX)
		    {
		      printf (" too many coloured pixels per proc\n");
		      printf (" the execution is terminated\n");
		      net->err = MPI_Abort (MPI_COMM_WORLD, 1);
		    }
		  if (rt->col_pixels == rt->col_pixels_max)
		    {
		      rt->col_pixels_max <<= 1;
		      rt->col_pixel_recv = (ColPixel *)realloc(rt->col_pixel_recv,
							       sizeof(ColPixel) * rt->col_pixels_max);
		    }
		  *col_pixel_id_p = rt->col_pixels;
		  
		  col_pixel_p = &rt->col_pixel_send[ *col_pixel_id_p ];
		  col_pixel_p->r = 0.F;
		  col_pixel_p->g = 0.F;
		  col_pixel_p->b = 0.F;
		  col_pixel_p->t = 1.e+30F;
		  col_pixel_p->i = (short int)i;
		  col_pixel_p->j = (short int)j;
		  ++rt->col_pixels;
		}
	      else
		{
		  col_pixel_p = &rt->col_pixel_send[ *col_pixel_id_p ];
		}
	      if (rt->is_isosurface)
		{
		  rt->t_min += t;
		  
		  if (rt->t_min < col_pixel_p->t)
		    {
		      col_pixel_p->r = ray.col_r;
		      col_pixel_p->g = ray.col_g;
		      col_pixel_p->b = ray.col_b;
		      col_pixel_p->t = rt->t_min;
		    }
		}
	      else
		{
		  col_pixel_p->r += ray.col_r;
		  col_pixel_p->g += ray.col_g;
		  col_pixel_p->b += ray.col_b;
		}
	    }
	  par3x += par1x;
	  par3y += par1y;
	  par3z += par1z;
	}
    }
  
  // here, intra-machine communications are handled through a binary
  // tree pattern and parallel pairwise syncronous communications. The
  // master processor of the current machine gets the sub-images of
  // all the processors of that machine. Inter-communications, needed
  // if the number of machines is greater than one, take place in the
  // routine rtRayTracingB.
  
  memcpy (rt->col_pixel_recv, rt->col_pixel_send, rt->col_pixels * sizeof(ColPixel));

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
	      net->err = MPI_Send (&rt->col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
	      
	      if (rt->col_pixels > 0)
		{
		  net->err = MPI_Send (&rt->col_pixel_send, rt->col_pixels, MPI_col_pixel_type, recv_id, 20, MPI_COMM_WORLD);
		}
	    }
	  else
	    {
	      net->err = MPI_Recv (&col_pixels, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD,
				   net->status);
	      
	      if (col_pixels > 0)
		{
		  net->err = MPI_Recv (&rt->col_pixel_send, col_pixels, MPI_col_pixel_type, send_id, 20, MPI_COMM_WORLD,
				       net->status);
		}
	      for (n = 0; n < col_pixels; n++)
		{
		  col_pixel1 = &rt->col_pixel_send[ n ];
		  i = col_pixel1->i;
		  j = col_pixel1->j;
		  
		  if (*(col_pixel_id_p = &rt->col_pixel_id[ i * pixels_y + j ]) == -1)
		    {
		      if (rt->col_pixels == rt->col_pixels_max)
			{
			  rt->col_pixels_max <<= 1;
			  rt->col_pixel_recv = (ColPixel *)realloc(rt->col_pixel_recv,
								   sizeof(ColPixel) * rt->col_pixels_max);
			}
		      col_pixel2 = &rt->col_pixel_recv[ *col_pixel_id_p = rt->col_pixels ];
		      
		      memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		      ++rt->col_pixels;
		    }
		  else
		    {
		      col_pixel2 = &rt->col_pixel_recv[ *col_pixel_id_p ];
		      
		      if (rt->is_isosurface)
			{
			  if (col_pixel1->t < col_pixel2->t)
			    {
			      col_pixel2->r = col_pixel1->r;
			      col_pixel2->g = col_pixel1->g;
			      col_pixel2->b = col_pixel1->b;
			      col_pixel2->t = col_pixel1->t;
			    }
			}
		      else
			{
			  col_pixel2->r += col_pixel1->r;
			  col_pixel2->g += col_pixel1->g;
			  col_pixel2->b += col_pixel1->b;
			}
		    }
		}
	    }
	  if (m < net->procs_per_machine[ machine_id ])
	    {
	      if (net->id == recv_id)
		{
		  memcpy (rt->col_pixel_send, rt->col_pixel_recv,
			  rt->col_pixels * sizeof(ColPixel));
		}
	    }
	  recv_id += comm_inc << 1;
	}
      comm_inc <<= 1;
    }
  
  if (net->machines == 1 || net->id != master_proc_id) return;
  
  // non-blocking inter-machine communications of sub-images begin here
  
  for (k = 0; k < 4 * pixels_x * pixels_y;)
    {
      pixel_color_to_send[ k   ] = 0.F;
      pixel_color_to_send[ k+1 ] = 0.F;
      pixel_color_to_send[ k+2 ] = 0.F;
      pixel_color_to_send[ k+3 ] = 1.e+30F;
      k += 4;
    }
  for (n = 0; n < rt->col_pixels; n++)
    {
      col_pixel_p = &rt->col_pixel_recv[ n ];
      
      i = col_pixel_p->i;
      j = col_pixel_p->j;
      
      k = 4 * (i * pixels_y + j);
      
      pixel_color_to_send[ k   ] += col_pixel_p->r;
      pixel_color_to_send[ k+1 ] += col_pixel_p->g;
      pixel_color_to_send[ k+2 ] += col_pixel_p->b;
      pixel_color_to_send[ k+3 ]  = col_pixel_p->t;
    }
  if (net->id != 0)
    {
      recv_id = 0;
      
      net->err = MPI_Issend (pixel_color_to_send,
			     4 * pixels_x * pixels_y, MPI_FLOAT,
			     recv_id, 30, MPI_COMM_WORLD,
			     &net->req[ 2 ][ net->id * net->procs + recv_id ]);
    }
  else
    {
      send_id = net->procs_per_machine[ net->id ];
      
      for (m = 1; m < net->machines; m++)
	{
	  net->err = MPI_Irecv (&pixel_color_to_recv[ (m-1) * (4 * SCREEN_SIZE_MAX) ],
				4 * pixels_x * pixels_y, MPI_FLOAT,
				send_id, 30, MPI_COMM_WORLD,
				&net->req[ 2 ][ (net->id + net->procs) * net->procs + send_id ]);
	  send_id += net->procs_per_machine[ m ];
	}
    }
}


void rtRayTracingB (void (*AbsorptionCoefficients) (float flow_field_data, float t1, float t2,
						    float cutoff, float *r, float *g, float *b),
		   Net *net, RT *rt)
{
  // here, the intra-machine communications take place and the RG
  // buffers or the output file are set
  
  float factor;
  float r, g, b;
  
  int pixels_x, pixels_y;
  int i, j, k;
  int m, n;
  int send_id, recv_id;
  int master_proc_id;
  int offset;
  
  short int pixel_i, pixel_j;
  
  unsigned char pixel_r, pixel_g, pixel_b;

  ColPixel *col_pixel_p;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  if (net->machines > 1)
    {
      master_proc_id = 0;
      
      for (m = 0; m < net->machine_id[ net->id ]; m++)
	{
	  master_proc_id += net->procs_per_machine[ m ];
	}
      if (net->id != master_proc_id) return;
      
      if (net->id != 0)
	{
	  recv_id = 0;
	  
	  net->err = MPI_Wait (&net->req[ 2 ][ net->id * net->procs + recv_id ], net->status);
	}
      else
	{
	  send_id = net->procs_per_machine[ net->id ];
	  
	  for (m = 1; m < net->machines; m++)
	    {
	      net->err = MPI_Wait (&net->req[ 2 ][ (net->id + net->procs) * net->procs + send_id ], net->status);
	      
	      offset = (m-1) * (4 * SCREEN_SIZE_MAX);
	      
	      for (k = 0; k < 4 * pixels_x * pixels_y;)
		{
		  if (rt->is_isosurface)
		    {
		      if (pixel_color_to_recv[ offset+3 ] < pixel_color_to_send[ k+3 ])
			{
			  pixel_color_to_send[ k   ] = pixel_color_to_recv[ offset   ];
			  pixel_color_to_send[ k+1 ] = pixel_color_to_recv[ offset+1 ];
			  pixel_color_to_send[ k+2 ] = pixel_color_to_recv[ offset+2 ];
			  pixel_color_to_send[ k+3 ] = pixel_color_to_recv[ offset+3 ];
			}
		    }
		  else
		    {
		      pixel_color_to_send[ k   ] += pixel_color_to_recv[ offset   ];
		      pixel_color_to_send[ k+1 ] += pixel_color_to_recv[ offset+1 ];
		      pixel_color_to_send[ k+2 ] += pixel_color_to_recv[ offset+2 ];
		    }
		  k += 4;
		  offset += 4;
		}
	      send_id += net->procs_per_machine[ m ];
	    }
	}
    }
  
  if (net->id != 0) return;
  
#ifndef BENCH
  if (rt->is_isosurface)
    {
      factor = 255.F;
    }
  else
    {
      factor = 255.F * rt->absorption_factor;
    }
#ifdef RG
  
  int bytes_per_pixel = sizeof(unsigned char) * 3;
  
  
  if (net->machines == 1)
    {
      for (n = 0; n < (pixels_x * pixels_y) * 3; n++)
	{
	  pixel_data[ n ] = 255;
	}
      for (n = 0; n < rt->col_pixels; n++)
	{
	  col_pixel_p = &rt->col_pixel_recv[ n ];
	  
	  r = col_pixel_p->r;
	  g = col_pixel_p->g;
	  b = col_pixel_p->b;
	  i = col_pixel_p->i;
	  j = col_pixel_p->j;
	  
	  k = (i * pixels_y + j) * 3;
	  
	  pixel_data[ k   ] = (unsigned char)max(0, min(255, (int)(255.F - factor * r)));
	  pixel_data[ k+1 ] = (unsigned char)max(0, min(255, (int)(255.F - factor * g)));
	  pixel_data[ k+2 ] = (unsigned char)max(0, min(255, (int)(255.F - factor * b)));
	}
    }
  else
    {
      for (k = 0; k < pixels_x * pixels_y; k++)
	{
	  r = pixel_color_to_send[ (k<<2)   ];
	  g = pixel_color_to_send[ (k<<2)+1 ];
	  b = pixel_color_to_send[ (k<<2)+2 ];
	  
	  pixel_data[ (k*3)   ] = (unsigned char)max(0, min(255, (int)(255.F - factor * r)));
	  pixel_data[ (k*3)+1 ] = (unsigned char)max(0, min(255, (int)(255.F - factor * g)));
	  pixel_data[ (k*3)+2 ] = (unsigned char)max(0, min(255, (int)(255.F - factor * b)));
	}
    }
  
  int eVizret = eViz_RLE_writeMemory (pixel_data, pixels_x, pixels_y,
				      bytes_per_pixel, &compressed_frame_size, compressed_data);
  
#else // RG
  
  FILE *image_file;
  XDR	xdr_image_file;
  
  
  image_file = fopen (rt->image_file_name, "w");
  xdrstdio_create (&xdr_image_file, image_file, XDR_ENCODE);
  
  xdr_int (&xdr_image_file, &pixels_x);
  xdr_int (&xdr_image_file, &pixels_y);
  
  if (net->machines == 1)
    {
      xdr_int (&xdr_image_file, &rt->col_pixels);
      
      for (n = 0; n < rt->col_pixels; n++)
	{
	  col_pixel_p = &rt->col_pixel_recv[ n ];
	  
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
    }
  else
    {
      int dummy = -1;
      
      xdr_int (&xdr_image_file, &dummy);
      
      for (k = 0; k < pixels_x * pixels_y; k++)
	{
	  r = pixel_color_to_send[ (k<<2)   ];
	  g = pixel_color_to_send[ (k<<2)+1 ];
	  b = pixel_color_to_send[ (k<<2)+2 ];
	  
	  pixel_r = (unsigned char)max(0, min(255, (int)(255.F - factor * r)));
	  pixel_g = (unsigned char)max(0, min(255, (int)(255.F - factor * g)));
	  pixel_b = (unsigned char)max(0, min(255, (int)(255.F - factor * b)));
	  
	  xdr_u_char (&xdr_image_file, &pixel_r);
	  xdr_u_char (&xdr_image_file, &pixel_g);
	  xdr_u_char (&xdr_image_file, &pixel_b);
	}
    }
  
  xdr_destroy (&xdr_image_file);
  fclose (image_file);
  
#endif // RG
#endif // BENCH
}


void rtRotate (float sin_1, float cos_1,
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


void rtProjection (float ortho_x, float ortho_y,
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
  
  rtRotate (viewpoint.sin_1, viewpoint.cos_1,
	    viewpoint.sin_2, viewpoint.cos_2,
	    screen.max_x, 0.0F, 0.0F,
	    &screen.dir1x, &screen.dir1y, &screen.dir1z);
  
  rtRotate (viewpoint.sin_1, viewpoint.cos_1,
	    viewpoint.sin_2, viewpoint.cos_2,
	    0.0F, screen.max_y, 0.0F,
	    &screen.dir2x, &screen.dir2y, &screen.dir2z);
  
  screen.par_x = (2.F * screen.max_x) / (float)pixels_x;
  screen.par_y = (2.F * screen.max_y) / (float)pixels_y;
}

#ifdef STEER
void rtReadParameters (char *parameters_file_name, Net *net, RT *rt, SteerParams *steer)
#else
void rtReadParameters (char *parameters_file_name, Net *net, RT *rt)
#endif
{
  FILE *parameters_file;
  
  float par_to_send[ 16 ];
  float ctr_x, ctr_y, ctr_z;
  float longitude, latitude;
  float zoom;
  float density_max, velocity_max, stress_max;
  
  int pixels_x, pixels_y;
  
  if (net->id == 0)
    {
      parameters_file = fopen (parameters_file_name, "r");
      
      fscanf (parameters_file, "%i \n", &pixels_x);
      fscanf (parameters_file, "%i \n", &pixels_y);
      fscanf (parameters_file, "%e \n", &ctr_x);
      fscanf (parameters_file, "%e \n", &ctr_y);
      fscanf (parameters_file, "%e \n", &ctr_z);
      fscanf (parameters_file, "%e \n", &longitude);
      fscanf (parameters_file, "%e \n", &latitude);
      fscanf (parameters_file, "%e \n", &zoom);
      
      fscanf (parameters_file, "%i \n", &rt->image_frequency);
      fscanf (parameters_file, "%i \n", &rt->flow_field_type);
      fscanf (parameters_file, "%i \n", &rt->is_isosurface);
      fscanf (parameters_file, "%e \n", &rt->absorption_factor);
      fscanf (parameters_file, "%e \n", &rt->cutoff);
      fscanf (parameters_file, "%e \n", &density_max);
      fscanf (parameters_file, "%e \n", &velocity_max);
      fscanf (parameters_file, "%e \n", &stress_max);
      
      fclose (parameters_file);
  
      par_to_send[  0 ] = 0.1 + (float)pixels_x;
      par_to_send[  1 ] = 0.1 + (float)pixels_y;
      par_to_send[  2 ] = ctr_x;
      par_to_send[  3 ] = ctr_y;
      par_to_send[  4 ] = ctr_z;
      par_to_send[  5 ] = longitude;
      par_to_send[  6 ] = latitude;
      par_to_send[  7 ] = zoom;
      par_to_send[  8 ] = 0.1 + (float)rt->image_frequency;
      par_to_send[  9 ] = 0.1 + (float)rt->flow_field_type;
      par_to_send[ 10 ] = 0.1 + (float)rt->is_isosurface;
      par_to_send[ 11 ] = rt->absorption_factor;
      par_to_send[ 12 ] = rt->cutoff;
      par_to_send[ 13 ] = density_max;
      par_to_send[ 14 ] = velocity_max;
      par_to_send[ 15 ] = stress_max;
    }
  net->err = MPI_Bcast (par_to_send, 16, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
  pixels_x              = (int)par_to_send[  0 ];
  pixels_y              = (int)par_to_send[  1 ];
  ctr_x                 =      par_to_send[  2 ];
  ctr_y                 =      par_to_send[  3 ];
  ctr_z                 =      par_to_send[  4 ];
  longitude             =      par_to_send[  5 ];
  latitude              =      par_to_send[  6 ];
  zoom                  =      par_to_send[  7 ];
  rt->image_frequency   = (int)par_to_send[  8 ];
  rt->flow_field_type   = (int)par_to_send[  9 ];
  rt->is_isosurface     = (int)par_to_send[ 10 ];
  rt->absorption_factor =      par_to_send[ 11 ];
  rt->cutoff            =      par_to_send[ 12 ];
  density_max           =      par_to_send[ 13 ];
  velocity_max          =      par_to_send[ 14 ];
  stress_max            =      par_to_send[ 15 ];
  
  rtProjection (0.5F * rt->system_size, 0.5F * rt->system_size,
		pixels_x, pixels_y,
		ctr_x, ctr_y, ctr_z,
		5.F * rt->system_size,
		longitude, latitude,
		0.5F * (5.F * rt->system_size),
		zoom);
  
  rt->flow_field_value_max_inv[ DENSITY  ] = 1.F / density_max;
  rt->flow_field_value_max_inv[ VELOCITY ] = 1.F / velocity_max;
  rt->flow_field_value_max_inv[ STRESS   ] = 1.F / stress_max;
  
#ifdef STEER
  // set up the ReG struct
  
  steer->pixels_x = pixels_x;
  steer->pixels_y = pixels_y;
  steer->longitude = longitude;
  steer->latitude = latitude;
  steer->zoom = zoom;
  steer->image_freq = rt->image_frequency;
  steer->flow_field_type = rt->flow_field_type;
  steer->is_isosurface = rt->is_isosurface;
  steer->abs_factor = rt->absorption_factor;
  steer->cutoff = rt->cutoff;
  steer->max_density = density_max;
  steer->max_velocity = velocity_max;
  steer->max_stress = stress_max;
#endif
}


#ifdef STEER
void rtUpdateParameters (RT *rt, SteerParams *steer)
{
  // update rt params
  
  rtProjection (0.5F * rt->system_size, 0.5F * rt->system_size,
		steer->pixels_x, steer->pixels_y,
		steer->ctr_x, steer->ctr_y, steer->ctr_z,
		5.F * rt->system_size,
		steer->longitude, steer->latitude,
		0.5F * (5.F * rt->system_size),
		steer->zoom);
  
  rt->image_frequency = steer->image_freq;
  rt->flow_field_type = steer->flow_field_type;
  rt->is_isosurface = steer->is_isosurface;
  rt->absorption_factor = steer->abs_factor;
  rt->cutoff = steer->cutoff;
  
  rt->flow_field_value_max_inv[ DENSITY  ] = 1.F / steer->max_density;
  rt->flow_field_value_max_inv[ VELOCITY ] = 1.F / steer->max_velocity;
  rt->flow_field_value_max_inv[ STRESS   ] = 1.F / steer->max_stress;
}
#endif


void rtInit (char *image_file_name, Net *net, RT *rt)
{
  int col_pixel_count = 6;
  int col_pixel_blocklengths[6] = {1, 1, 1, 1, 1, 1};
  MPI_Datatype col_pixel_types[6] = {MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL,
				     MPI_SHORT, MPI_SHORT};

  // calculate displacements
  MPI_Aint col_pixel_disps[6];
  
  
  rt->image_file_name = image_file_name;
  
  rt->col_pixels_max = 512 * 512;

  //rt->col_pixel_send = (ColPixel *)malloc(sizeof(ColPixel) * rt->col_pixels_max);

  rt->col_pixel_recv = (ColPixel *)malloc(sizeof(ColPixel) * rt->col_pixels_max);
  
  rt->pixels_max = 512 * 512;
  rt->col_pixel_id = (int *)malloc(sizeof(int) * rt->pixels_max);
  
  rt->col_pixels_recv = (int *)malloc(sizeof(int) * net->procs);
  
#ifdef RG
  if (net->id == 0)
    {
      pixel_data = (unsigned char *)malloc(sizeof(unsigned char) * 3 * rt->pixels_max);
      compressed_data = (unsigned char *)malloc(sizeof(unsigned char) * 3 * rt->pixels_max);
    }
#endif
  
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

void rtEnd (Net *net, RT *rt)
{
  free(rt->col_pixels_recv);
  //free(rt->col_pixel_send);
  free(rt->col_pixel_id);
  free(rt->col_pixel_recv);
  free(rt->cluster);
  
#ifdef RG
  if (net->id == 0)
    {
      free(compressed_data);
      free(pixel_data);
    }
#endif // RG
}
