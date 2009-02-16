#include "segmentation.h"


void segFromVoxelCoordsToSiteCoords (short int voxel[3], short int site[3], Vis *vis)
{
  for (int l = 0; l < 3; l++)
    site[l] = (int)(voxel[l] * vis->scale[l]);
}


void segFromSiteCoordsToVoxelCoords (short int site[3], short int voxel[3], Vis *vis)
{
  for (int l = 0; l < 3; l++)
    voxel[l] = (int)(site[l] * vis->inv_scale[l]);
}


void segFromSiteCoordsToGridCoords (short int site[3], int *b_id, int *s_id, Vis *vis)
{
  int b[3], s[3];
  
  
  for (int l = 0; l < 3; l++)
    b[l] = site[l] >> SHIFT;
  
  *b_id = BlockId(b,vis->blocks);
  
  for (int l = 0; l < 3; l++)
    s[l] = site[l] - (b[l] << SHIFT);
  
  *s_id = (((s[0] << SHIFT) + s[1]) << SHIFT) + s[2];
}


Site *segSitePtr (short int site[3], Vis *vis)
{
  int b[3];
  
  for (int l = 0; l < 3; l++)
    b[l] = site[l] >> SHIFT;
  
  Block *block_p = &vis->block[ BlockId(b,vis->blocks) ];
  
  if (block_p->site == NULL)
    {
      return NULL;
    }
  else
    {
      return &block_p->site[ SiteId(site,b) ];
    }
}


void segSetSiteData (Site *site_p, unsigned int config, unsigned int label)
{
  site_p->cfg   = config;
  site_p->label = label;
}


Site *segStackSitePtr (Vis *vis)
{
  vis->stack_sites += SITES_PER_BLOCK;
  
  if (vis->stack_sites > vis->stack_sites_max)
    {
      vis->stack_sites -= SITES_PER_BLOCK;
      return NULL;
    }
  for (int m = 0; m < SITES_PER_BLOCK; m++)
    vis->stack_site[ vis->stack_sites - SITES_PER_BLOCK + m ].cfg = SOLID_TYPE;
  
  return &vis->stack_site[ vis->stack_sites - SITES_PER_BLOCK ];
}


int segIsSuperficialSite (short int site[3], Vis *vis)
{
  short int neigh[3];
  
  
  for (int dir = 0; dir < 14; dir++)
    {
      for (int l = 0; l < 3; l++)
	{
	  neigh[l] = site[l] + e[ dir*3+l ];
	  
	  if (neigh[l] == -1 || neigh[l] == vis->sites[l])
	    {
	      return SUCCESS;
	    }
	}
      Site *neigh_site_p = segSitePtr (neigh, vis);
      
      if (neigh_site_p == NULL || neigh_site_p->cfg == SOLID_TYPE)
	{
	  return SUCCESS;
	}
    }
  return !SUCCESS;
}


int segRayVsDisc (Triangle *t_p, float x1[3], float x2[3], float t_max, float *t)
{
  float px1[3], px2[3];
  float px3[3], px4[3];
  float x[2];
  
  
  for (int l = 0; l < 3; l++)
    {
      px1[l] = x1[l] - t_p->pos[l];
      px2[l] = x2[l] - t_p->pos[l];
    }
  AntiRotate (px1,
	      t_p->d.sin_longitude, t_p->d.cos_longitude,
	      t_p->d.sin_latitude,  t_p->d.cos_latitude,
	      px3);
  
  AntiRotate (px2,
	      t_p->d.sin_longitude, t_p->d.cos_longitude,
	      t_p->d.sin_latitude,  t_p->d.cos_latitude,
	      px4);
  
  if ((px3[2] > 0.0F && px4[2] > 0.0F) ||
      (px3[2] < 0.0F && px4[2] < 0.0F))
    {
      return !SUCCESS;
    }
  *t = px3[2] / (px3[2] - px4[2]);
  
  if (*t >= t_max)
    {
      return !SUCCESS;
    }
  x[0] = (1.0F - *t) * px3[0] + *t * px4[0];
  x[1] = (1.0F - *t) * px3[1] + *t * px4[1];
  
  if (x[0]*x[0] + x[1]*x[1] > t_p->d.r2)
    {
      return !SUCCESS;
    }
  return SUCCESS;
}


float segInterpolatedGrey (short int site[3], Vis *vis)
{
  if (vis->res_factor == 1)
    {
      if (site[0] >= vis->input_voxels[0] ||
	  site[1] >= vis->input_voxels[1] ||
	  site[2] >= vis->input_voxels[2])
	{
	  return 1.0e+30F;
	}
      return (float)vis->voxel[ VoxelId(site,vis->input_voxels) ];
    }
  float grey[2][2][2];
  float x[3];
  
  int voxel[2][3], v[3];
  int i, j, k, l;
  
  
  for (l = 0; l < 3; l++)
    {
      voxel[0][l] = (int)(x[l] = (float)site[l] * vis->inv_scale[l]);
      x[l] -= (float)voxel[0][l];
      voxel[1][l] = min(voxel[0][l] + 1, vis->input_voxels[l] - 1);
    }
  for (k = 0; k < 2; k++)
    {
      v[2] = voxel[k][2];
      
      for (j = 0; j < 2; j++)
        {
	  v[1] = voxel[j][1];
	  
          for (i = 0; i < 2; i++)
            {
	      v[0] = voxel[i][0];
	      
	      grey[k][j][i] = (float)vis->voxel[ VoxelId(v,vis->input_voxels) ];
	    }
	}
    }
  float grey_00x = (1.0F - x[0]) * grey[0][0][0] + x[0] * grey[0][0][1];
  float grey_01x = (1.0F - x[0]) * grey[0][1][0] + x[0] * grey[0][1][1];
  float grey_10x = (1.0F - x[0]) * grey[1][0][0] + x[0] * grey[1][0][1];
  float grey_11x = (1.0F - x[0]) * grey[1][1][0] + x[0] * grey[1][1][1];
  
  float grey_0y = (1.0F - x[1]) * grey_00x + x[1] * grey_01x;
  float grey_1y = (1.0F - x[1]) * grey_10x + x[1] * grey_11x;
  
  return (1.0F - x[2]) * grey_0y + x[2] * grey_1y;
}


int segEstimateNormal (short int site[3], float nor[], Vis *vis)
{
  float org[3];
  float x[3];
  float temp;
  
  int block_id, site_id;
  int is_front_advancing;
  int iters;
  int i, l;
  
  short int neigh[3];
  
  Site *site_p;
  
  Block *block_p;
  
  
  segFromSiteCoordsToGridCoords (site, &block_id, &site_id, vis);
  
  block_p = &vis->block[ block_id ];
  
  if (block_p->site == NULL)
    return !SUCCESS;
  
  for (i = 0; i < SITES_PER_BLOCK; i++)
    block_p->site[i].label = 0;
  
  block_p->site[ site_id ].label = 1;
  
  for (l = 0; l < 3; l++)
    {
      org[l] = (float)site[l];
      nor[l] = 0.F;
    }
  for (l = 0; l < 3; l++)
    vis->coord[A][0].x[l] = site[l];
  vis->coords[A] = 1;
  
  is_front_advancing = 1;
  iters = 0;
  
  while (++iters <= ITERS_MAX && is_front_advancing)
    {
      is_front_advancing = 0;
      vis->coords[B] = 0;
      
      for (i = 0; i < vis->coords[A]; i++)
	{
	  for (int dir = 0; dir < 14; dir++)
	    {
	      for (l = 0; l < 3; l++)
		{
		  neigh[l] = vis->coord[A][i].x[l] + e[ dir*3+l ];
		  
		  if (neigh[l] == -1 || neigh[l] == vis->sites[l])
		    {
		      l = 10;
		      break;
		    }
		}
	      if (l >= 10) continue;
	      
	      segFromSiteCoordsToGridCoords (neigh, &block_id, &site_id, vis);
	      
	      if ((block_p = &vis->block[ block_id ])->site == NULL)
		{
		  continue;
		}
	      site_p = &block_p->site[ site_id ];
	      
	      if ((site_p->cfg & SITE_TYPE_MASK) != FLUID_TYPE ||
		  site_p->label == 1 ||
		  !segIsSuperficialSite (neigh, vis))
		{
		  continue;
		}
	      if (vis->coords[B] >= COORD_BUFFER_SIZE_B)
		{
		  return !SUCCESS;
		}
	      is_front_advancing = 1;
	      
	      site_p->label = 1;
	      
	      for (l = 0; l < 3; l++)
		vis->coord[B][ vis->coords[B] ].x[l] = neigh[l];
	      ++vis->coords[B];
	      
	      for (l = 0; l < 3; l++)
		x[l] = (float)neigh[l] - org[l];
	      
	      temp = 1.0 / sqrtf(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	      
	      for (l = 0; l < 3; l++)
		nor[l] += x[l] * temp;
	    }
	}
      Coord *temp = vis->coord[A];
      vis->coord[A] = vis->coord[B];
      vis->coord[B] = temp;
      vis->coords[A] = vis->coords[B];
    }
  temp = 1.0 / sqrtf(nor[0]*nor[0] + nor[1]*nor[1] + nor[2]*nor[2]);
  
  for (l = 0; l < 3; l++)
    nor[l] *= temp;
  
  return SUCCESS;
}


void segEstimateDiameter (short int site[3], float nor[3], float *diameter, Vis *vis)
{
  float segment_length = 0.25F;
  float org[3], x[3];
  
  int iters;
  int l;
  
  short int neigh[3];
  
  Site *neigh_site_p;
  
  
  for (l = 0; l < 3; l++)
    {
      org[l] = (float)site[l];
      neigh[l] = (int)(org[l] + 2.0F * nor[l]);
    }
  if (neigh[0] < 0 || neigh[0] >= vis->sites[0] ||
      neigh[1] < 0 || neigh[1] >= vis->sites[1] ||
      neigh[2] < 0 || neigh[2] >= vis->sites[2])
    {
      for (l = 0; l < 3; l++)
	nor[l] = -nor[l];
    }
  else
    {
      neigh_site_p = segSitePtr (neigh, vis);
      
      if (neigh_site_p == NULL || neigh_site_p->cfg == SOLID_TYPE)
	{
	  for (l = 0; l < 3; l++)
	    nor[l] = -nor[l];
	}
    }
  for (l = 0; l < 3; l++)
    {
      nor[l] *= segment_length;
      x[l] = org[l];
    }
  iters = 0;
  
  neigh_site_p = segSitePtr (site, vis);
  
  while (neigh_site_p != NULL && neigh_site_p->cfg != SOLID_TYPE)
    {
      ++iters;
      
      for (l = 0; l < 3; l++)
	{
	  x[l] += nor[l];
	  neigh[l] = (int)x[l];
	}
      if (neigh[0] < 0 || neigh[0] >= vis->sites[0] ||
	  neigh[1] < 0 || neigh[1] >= vis->sites[1] ||
	  neigh[2] < 0 || neigh[2] >= vis->sites[2])
	{
	  neigh_site_p = NULL;
	}
      else
	{
	  neigh_site_p = segSitePtr (neigh, vis);
	}
    }
  --iters;
  
  if (iters == 0)
    {
      *diameter = -1.0;
    }
  else
    {
      *diameter = iters * segment_length;
    }
}


int segCreateOptimisedTriangle (unsigned int site_type, short int site[3], Vis *vis)
{
  float triangle_factor = 2.F;
  float org[3], nor[3];
  float x[3];
  float longitude, latitude;
  float triangle_size;
  float diameter;
  
  int k, l, m, n;
  int iters, iters_max;
  int longitude_id, latitude_id;
  int t_id;
  
  short int neigh[3];
  
  Site *site_p, *neigh_site_p;
  
  Triangle *t_p;
  
  
  if (vis->boundary[ site_type ].triangles == (int)(1U << BOUNDARY_ID_BITS))
    {
      return -1;
    }
  if (segEstimateNormal (site, nor, vis) == !SUCCESS)
    {
      return -1;
    }
  segEstimateDiameter (site, nor, &diameter, vis);
  
  if (diameter < 0.F)
    {
      return -1;
    }
  for (l = 0; l < 3; l++)
    {
      org[l] = (float)site[l];
      org[l] = 0.5F * (org[l] + (org[l] + diameter * nor[l]));
      site[l] = (int)org[l];
    }
  site_p = segSitePtr (site, vis);
  
  iters_max = 0;
  longitude_id = 0;
  latitude_id = 0;
  
  latitude = 0.F;
  
  for (n = 0; n < 180; n++)
    {
      longitude = 0.F;
      
      for (m = 0; m < 360; m++)
	{
	  x[0] = 0.0F;
	  x[1] = 0.0F;
	  x[2] = 1.0F;
	  
	  Rotate (x, longitude, latitude, nor);
	  
	  for (l = 0; l < 3; l++)
	    x[l] = org[l];
	  
	  iters = 0;
	  neigh_site_p = site_p;
	  
	  while (neigh_site_p != NULL && neigh_site_p->cfg != SOLID_TYPE)
	    {
	      ++iters;
	      
	      for (l = 0; l < 3; l++)
		{
		  x[l] += nor[l];
		  neigh[l] = (int)x[l];
		}
	      if (neigh[0] < 0 || neigh[0] >= vis->sites[0] ||
		  neigh[1] < 0 || neigh[1] >= vis->sites[1] ||
		  neigh[2] < 0 || neigh[2] >= vis->sites[2])
		{
		  neigh_site_p = NULL;
		}
	      else
		{
		  neigh_site_p = segSitePtr (neigh, vis);
		}
	    }
	  --iters;
	  
	  if (iters > iters_max)
	    {
	      iters_max = iters;
	      longitude_id = m;
	      latitude_id = n;
	    }
	  longitude += DEG_TO_RAD;
	}
      latitude += 2.0F * DEG_TO_RAD;
    }
  if (iters_max == 0)
    {
      return -1;
    }
  longitude = longitude_id * DEG_TO_RAD;
  latitude = latitude_id * 2.0F * DEG_TO_RAD;
  
  t_id = vis->boundary[ site_type ].triangles;
  t_p = &vis->boundary[ site_type ].triangle[ t_id ];

  triangle_size = triangle_factor * diameter;
  
  for (l = 0; l < 3; l++)
    org[l] -= vis->half_dim[l];
  
  x[0] = 0.0F;
  x[1] = 2.0F * triangle_size;
  x[2] = 0.0F;
  
  Rotate (x, longitude, latitude, t_p->v[0].pos);
  
  x[0] = -(triangle_size * sqrtf(3.0F));
  x[1] = -triangle_size;
  x[2] = 0.0F;
  
  Rotate (x, longitude, latitude, t_p->v[1].pos);
  
  x[0] = +(triangle_size / sqrtf(3.0F));
  x[1] = -triangle_size;
  x[2] = 0.0F;
  
  Rotate (x, longitude, latitude, t_p->v[2].pos);
  
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++)
      t_p->v[k].pos[l] += org[l];
  
  t_p->pressure_avg = 80.0F;
  t_p->pressure_amp = 0.0F;
  t_p->pressure_phs = 0.0F;
  
  t_p->normal_sign = 1;
  
  editCalculateTriangleData (t_p);
  
  ++vis->boundary[ site_type ].triangles;
  
  return t_id;
}


int segIsSegmentIntercepted (short int site[3], int dir, Site *site_p, Vis *vis)
{
  float x1[3], x2[3];
  float t, t_max;
  
  int b_id, t_id;
  
  
  for (int l = 0; l < 3; l++)
    x2[l] = (x1[l] = site[l] - vis->half_dim[l]) + e[ dir*3+l ];
  
  t_max = 1.0;
  b_id = -1;
  
  for (int n = 0; n < BOUNDARIES; n++)
    for (int m = 0; m < vis->boundary[ n ].triangles; m++)
      {
	if (segRayVsDisc (&vis->boundary[n].triangle[m], x1, x2, t_max, &t))
	{
	  t_max = t;
	  b_id = n;
	  t_id = m;
	}
    }
  if (b_id == -1)
    {
      return !SUCCESS;
    }
  else
    {
      if (b_id == INLET_BOUNDARY)
	{
	  site_p->cfg = INLET_TYPE | (t_id << BOUNDARY_ID_SHIFT);
	}
      else if (b_id == OUTLET_BOUNDARY)
	{
	  site_p->cfg = OUTLET_TYPE | (t_id << BOUNDARY_ID_SHIFT);
	}
      return SUCCESS;
    }
}


int segSegmentation (Vis *vis)
{
  float seconds = myClock ();
  
  int block_id, site_id;
  int is_front_advancing, iters;
  int i, l;
  
  short int site[3], neigh[3];
  
  Site *site_p;
  
  Block *block_p;
  
  
  for (i = 0; i < vis->tot_blocks; i++)
    {
      vis->block[i].site = NULL;
    }
  vis->tot_sites = 0;
  vis->stack_sites = 0;
  vis->coords[A] = 0;
  vis->coords[C] = 0;
  
  segFromVoxelCoordsToSiteCoords (vis->selected_voxel, site, vis);
  
  if (segInterpolatedGrey (site, vis) < vis->selected_grey)
    {
      vis->segmentation_time = myClock () - seconds;
      return SUCCESS;
    }
  segFromSiteCoordsToGridCoords (site, &block_id, &site_id, vis);
  
  (block_p = &vis->block[ block_id ])->site = segStackSitePtr (vis);
  
  segSetSiteData (&block_p->site[ site_id ], FLUID_TYPE, 0);
  
  vis->tot_sites = 1;
  iters = 0;
  
  for (l = 0; l < 3; l++)
    vis->coord[A][0].x[l] = site[l];
  vis->coord[A][0].iters = iters;
  vis->coords[A] = 1;
  
  is_front_advancing = 1;
  
  while (is_front_advancing)
    {
      is_front_advancing = 0;
      ++iters;
      vis->coords[B] = 0;
      
      for (i = 0; i < vis->coords[A]; i++)
	{
	  for (l = 0; l < 3; l++)
	    site[l] = vis->coord[A][i].x[l];
	  
	  site_p = segSitePtr (site, vis);
	  
	  for (int dir = 0; dir < 14; dir++)
	    {
	      for (l = 0; l < 3; l++)
		{
		  neigh[l] = site[l] + e[ dir*3+l ];
		  
		  if (neigh[l] == -1 || neigh[l] == vis->sites[l])
		    {
		      l = 10;
		      break;
		    }
		}
	      if (l >= 10) continue;
	      
	      segFromSiteCoordsToGridCoords (neigh, &block_id, &site_id, vis);
	      
	      block_p = &vis->block[ block_id ];
	      
	      if (block_p->site != NULL && block_p->site[ site_id ].cfg != SOLID_TYPE)
		{
		  continue;
		}
	      if (segInterpolatedGrey (neigh, vis) < vis->selected_grey ||
		  segIsSegmentIntercepted (site, dir, site_p, vis) == SUCCESS)
		{
		  continue;
		}
	      if (vis->coords[B] >= COORD_BUFFER_SIZE_B) return !SUCCESS;
	      
	      is_front_advancing = 1;
	      
	      if (block_p->site == NULL)
		{
		  if ((block_p->site = segStackSitePtr (vis)) == NULL)
		    {
		      return !SUCCESS;
		    }
		}
	      segSetSiteData (&block_p->site[ site_id ], FLUID_TYPE, 0);
	      
	      ++vis->tot_sites;
	      
	      for (l = 0; l < 3; l++)
		vis->coord[B][ vis->coords[B] ].x[l] = neigh[l];
	      vis->coord[B][ vis->coords[B] ].iters = iters;
	      ++vis->coords[B];
	    }
	}
      for (i = 0; i < vis->coords[A]; i++)
	{
	  if (segIsSuperficialSite (vis->coord[A][i].x, vis))
	    {
	      if (vis->coords[C] >= COORD_BUFFER_SIZE_C) return !SUCCESS;
	      
	      memcpy (&vis->coord[C][ vis->coords[C] ], &vis->coord[A][i], sizeof(Coord));
	      ++vis->coords[C];
	    }
	}
      Coord *temp = vis->coord[A];
      vis->coord[A] = vis->coord[B];
      vis->coord[B] = temp;
      vis->coords[A] = vis->coords[B];
    }
  vis->segmentation_time = myClock () - seconds;
  
  return SUCCESS;
}


int segUpdateSegmentation (Vis *vis)
{
  float seconds = myClock ();
  
  int block_id, site_id;
  int is_front_advancing;
  int i, l;
  
  short int site[3], neigh[3];
  
  Site *site_p;
  
  Block *block_p;
  
  
  memcpy (&vis->coord[A][0], &vis->coord[C][0], vis->coords[C] * sizeof(Coord));
  vis->coords[A] = vis->coords[C];
  vis->coords[C] = 0;
  
  is_front_advancing = 1;
  
  while (is_front_advancing)
    {
      is_front_advancing = 0;
      vis->coords[B] = 0;
      
      for (i = 0; i < vis->coords[A]; i++)
	{
	  for (l = 0; l < 3; l++)
	    site[l] = vis->coord[A][i].x[l];
	  
	  site_p = segSitePtr (site, vis);
	  
	  for (int dir = 0; dir < 14; dir++)
	    {
	      for (l = 0; l < 3; l++)
		{
		  neigh[l] = site[l] + e[ dir*3+l ];
		  
		  if (neigh[l] == -1 || neigh[l] == vis->sites[l])
		    {
		      l = 10;
		      break;
		    }
		}
	      if (l >= 10) continue;
	      
	      segFromSiteCoordsToGridCoords (neigh, &block_id, &site_id, vis);
	      
	      block_p = &vis->block[ block_id ];
	      
	      if (block_p->site != NULL && block_p->site[ site_id ].cfg != SOLID_TYPE)
		{
		  continue;
		}
	      if (segInterpolatedGrey (neigh, vis) < vis->selected_grey ||
		  segIsSegmentIntercepted (site, dir, site_p, vis) == SUCCESS)
		{
		  continue;
		}
	      if (vis->coords[B] >= COORD_BUFFER_SIZE_B) return !SUCCESS;
	      
	      is_front_advancing = 1;
	      
	      if (block_p->site == NULL)
		{
		  if ((block_p->site = segStackSitePtr (vis)) == NULL)
		    {
		      return !SUCCESS;
		    }
		}
	      segSetSiteData (&block_p->site[ site_id ], FLUID_TYPE, 0);
	      
	      ++vis->tot_sites;
	      
	      for (l = 0; l < 3; l++)
		vis->coord[B][ vis->coords[B] ].x[l] = neigh[l];
	      vis->coord[B][ vis->coords[B] ].iters = vis->coord[A][i].iters;
	      ++vis->coords[B];
	    }
	}
      for (i = 0; i < vis->coords[A]; i++)
	{
	  if (segIsSuperficialSite (vis->coord[A][i].x, vis))
	    {
	      if (vis->coords[C] >= COORD_BUFFER_SIZE_C) return !SUCCESS;
	      
	      memcpy (&vis->coord[C][ vis->coords[C] ], &vis->coord[A][i], sizeof(Coord));
	      ++vis->coords[C];
	    }
	}
      Coord *temp = vis->coord[A];
      vis->coord[A] = vis->coord[B];
      vis->coord[B] = temp;
      vis->coords[A] = vis->coords[B];
    }
  vis->segmentation_time = myClock () - seconds;
  
  return SUCCESS;
}


void segSetBoundaryConfigurations (Vis *vis)
{
  int unknowns;
  int i, j, k;
  int l, m, n;
  
  unsigned int boundary_config, dir;
  
  short int site[3], neigh[3];
  
  Site *site_p, *neigh_site_p;
  
  Block *block_p;
  
  
  n = -1;
  
  for (i = 0; i < vis->blocks[0]; i++)
    for (j = 0; j < vis->blocks[1]; j++)
      for (k = 0; k < vis->blocks[2]; k++)
	{
	  if ((block_p = &vis->block[++n])->site == NULL) continue;
	  
	  m = -1;
	  
	  for (site[0] = i * BLOCK_SIZE; site[0] < (i+1) * BLOCK_SIZE; site[0]++)
	    for (site[1] = j * BLOCK_SIZE; site[1] < (j+1) * BLOCK_SIZE; site[1]++)
	      for (site[2] = k * BLOCK_SIZE; site[2] < (k+1) * BLOCK_SIZE; site[2]++)
		{
		  if ((site_p = &block_p->site[++m])->cfg == SOLID_TYPE) continue;
		  
		  boundary_config = 0U;
		  unknowns = 0;
		  
		  for (dir = 0U; dir < 14U; dir++)
		    {
		      for (l = 0; l < 3; l++)
			{
			  neigh[l] = site[l] + e[ dir*3+l ];
			  
			  if (neigh[l] == -1 || neigh[l] == vis->sites[l])
			    {
			      l = 10;
			      break;
			    }
			}
		      if (l >= 10) continue;
		      
		      neigh_site_p = segSitePtr (neigh, vis);
		      
		      if (!(neigh_site_p == NULL || neigh_site_p->cfg == SOLID_TYPE))
			{
			  boundary_config |= (1U << dir);
			}
		      else
			{
			  ++unknowns;
			}
		    }
		  if ((site_p->cfg & SITE_TYPE_MASK) == FLUID_TYPE)
		    {
		      if (unknowns != 0)
			{
			  site_p->cfg |= (boundary_config << BOUNDARY_CONFIG_SHIFT);
			}
		    }
		  else
		    {
		      site_p->cfg |= (boundary_config << BOUNDARY_CONFIG_SHIFT);
		    }
		}
	}
}
