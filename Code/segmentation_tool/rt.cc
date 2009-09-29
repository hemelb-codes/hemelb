#ifdef MESH
#include "rt.h"


int rtTriangleVsRay (MeshTriangle *triangle, Ray *ray, Hit *hit)
{
  double e1[3], e2[3];
  double x1[3], x2[3];
  double det;
  double v, w, t;
  
  int l;
  
  
  for (l = 0; l < 3; l++)
    e1[l] = triangle->v[1].pos[l] - triangle->v[0].pos[l];
  
  for (l = 0; l < 3; l++)
    e2[l] = triangle->v[2].pos[l] - triangle->v[0].pos[l];
  
  x1[0] = ray->dir[1] * e2[2] - ray->dir[2] * e2[1];
  x1[1] = ray->dir[2] * e2[0] - ray->dir[0] * e2[2];
  x1[2] = ray->dir[0] * e2[1] - ray->dir[1] * e2[0];
  
  det = ScalarProd (e1, x1);
  
  if (det > -EPSILON && det < EPSILON)
    {
      return !SUCCESS;
    }
  for (l = 0; l < 3; l++)
    x2[l] = ray->org[l] - triangle->v[0].pos[l];
  
  v = ScalarProd (x1, x2) / det;
  
  if (v < 0.0 || v > 1.0)
    {
      return !SUCCESS;
    }
  x1[0] = x2[1] * e1[2] - x2[2] * e1[1];
  x1[1] = x2[2] * e1[0] - x2[0] * e1[2];
  x1[2] = x2[0] * e1[1] - x2[1] * e1[0];
  
  w = ScalarProd (ray->dir, x1) / det;
  
  if (w < 0.0 || v + w > 1.0)
    {
      return !SUCCESS;
    }
  t = ScalarProd (e2, x1) / det;
  
  if (t < ray->t_near || t > hit->t)
    {
      return !SUCCESS;
    }
  for (l = 0; l < 3; l++)
    hit->pos[l] = ray->org[l] + t * ray->dir[l];
  
  hit->t = t;
  
  return SUCCESS;
}


void rtAABBvsRay (double aabb[], Ray *ray)
{
  double tx0, ty0, tz0;
  double tx1, ty1, tz1;
  
  
  if (ray->dir[0] < 0.0)
    {
      tx0 = (aabb[0] - ray->org[0]) / ray->dir[0];
      tx1 = (aabb[1] - ray->org[0]) / ray->dir[0];
    }
  else
    {
      tx0 = (aabb[1] - ray->org[0]) / ray->dir[0];
      tx1 = (aabb[0] - ray->org[0]) / ray->dir[0];
    }
  if (ray->dir[1] < 0.0)
    {
      ty0 = (aabb[2] - ray->org[1]) / ray->dir[1];
      ty1 = (aabb[3] - ray->org[1]) / ray->dir[1];
    }
  else
    {
      ty0 = (aabb[3] - ray->org[1]) / ray->dir[1];
      ty1 = (aabb[2] - ray->org[1]) / ray->dir[1];
    }
  if (ray->dir[2] < 0.0)
    {
      tz0 = (aabb[4] - ray->org[2]) / ray->dir[2];
      tz1 = (aabb[5] - ray->org[2]) / ray->dir[2];
    }
  else
    {
      tz0 = (aabb[5] - ray->org[2]) / ray->dir[2];
      tz1 = (aabb[4] - ray->org[2]) / ray->dir[2];
    }
  ray->t_near = fmax(tx0, fmax(ty0, tz0));
  ray->t_far  = fmin(tx1, fmin(ty1, tz1));
}


void rtInitRayTracing (Mesh *mesh)
{
  int aabb[3][2];
  int i[3], voxel_id;
  int l, m, n;
  
  
  for (n = 0; n < mesh->voxels[3]; n++)
    {
      mesh->voxel[n].triangle_id = NULL;
      mesh->voxel[n].triangles = 0;
    }
  for (n = 0; n < mesh->triangles; n++)
    {
      for (l = 0; l < 3; l++)
	{
	  aabb[l][0] =  1000000000;
	  aabb[l][1] = -1000000000;
	}
      for (m = 0; m < 3; m++)
	for (l = 0; l < 3; l++)
	  {
	    voxel_id = (int)((mesh->triangle[n].v[m].pos[l] + mesh->half_dim[l]) / mesh->voxel_size);
	    
	    aabb[l][0] = min(aabb[l][0], voxel_id);
	    aabb[l][1] = max(aabb[l][1], voxel_id);
	  }
      for (i[0] = aabb[0][0]; i[0] <= aabb[0][1]; i[0]++)
	for (i[1] = aabb[1][0]; i[1] <= aabb[1][1]; i[1]++)
	  for (i[2] = aabb[2][0]; i[2] <= aabb[2][1]; i[2]++)
	    ++mesh->voxel[ (i[0]*mesh->voxels[1]+i[1])*mesh->voxels[2]+i[2] ].triangles;
    }
  for (n = 0; n < mesh->voxels[3]; n++)
    {
      if (mesh->voxel[n].triangles == 0) continue;
      
      mesh->voxel[n].triangle_id = (int *)malloc(sizeof(int) * mesh->voxel[n].triangles);
      mesh->voxel[n].triangles = 0;
    }
  for (n = 0; n < mesh->triangles; n++)
    {
      for (l = 0; l < 3; l++)
	{
	  aabb[l][0] =  1000000000;
	  aabb[l][1] = -1000000000;
	}
      for (m = 0; m < 3; m++)
	for (l = 0; l < 3; l++)
	  {
	    voxel_id = (int)((mesh->triangle[n].v[m].pos[l] + mesh->half_dim[l]) / mesh->voxel_size);
	    
	    aabb[l][0] = min(aabb[l][0], voxel_id);
	    aabb[l][1] = max(aabb[l][1], voxel_id);
	  }
      for (i[0] = aabb[0][0]; i[0] <= aabb[0][1]; i[0]++)
	for (i[1] = aabb[1][0]; i[1] <= aabb[1][1]; i[1]++)
	  for (i[2] = aabb[2][0]; i[2] <= aabb[2][1]; i[2]++)
	    {
	      voxel_id = (i[0]*mesh->voxels[1]+i[1])*mesh->voxels[2]+i[2];
	      mesh->voxel[ voxel_id ].triangle_id[ mesh->voxel[voxel_id].triangles++ ] = n;
	    }
    }
}


int rtVoxelVsRay (Voxel *voxel, Ray *ray, Hit *hit, Mesh *mesh)
{
  hit->t = fmin(ray->t_far, ray->t_max);
  hit->triangle_id = -1;
  
  for (int m = 0; m < voxel->triangles; m++)
    {
      if (voxel->triangle_id[m] == hit->previous_triangle_id) continue;
      
      if (rtTriangleVsRay (&mesh->triangle[ voxel->triangle_id[m] ], ray, hit))
	{
	  hit->triangle_id = voxel->triangle_id[m];
	}
    }
  if (hit->triangle_id != -1)
    {
      return SUCCESS;
    }
  else
    {
      return !SUCCESS;
    }
}


int rtTraceRay (Ray *ray, Hit *hit, Mesh *mesh)
{
  double t_max[3], t_delta[3];
  
  int i[3], inc[3];
  int l;
  
  Voxel *voxel;
  
  
  for (l = 0; l < 3; l++)
    i[l] = max(0, min(mesh->voxels[l] - 1,
		      (int)((ray->org[l] + ray->t_near * ray->dir[l] + mesh->half_dim[l]) / mesh->voxel_size)));
 
  if (ray->dir[0] < 0.0)
    {
      inc[0] = -1;
      t_max[0] = (i[0] * mesh->voxel_size - mesh->half_dim[0] - ray->org[0]) / ray->dir[0];
    }
  else
    {
      inc[0] = 1;
      t_max[0] = ((i[0] + 1) * mesh->voxel_size - mesh->half_dim[0] - ray->org[0]) / ray->dir[0];
    }
  if (ray->dir[1] < 0.0)
    {
      inc[1] = -1;
      t_max[1] = (i[1] * mesh->voxel_size - mesh->half_dim[1] - ray->org[1]) / ray->dir[1];
    }
  else
    {
      inc[1] = 1;
      t_max[1] = ((i[1] + 1) * mesh->voxel_size - mesh->half_dim[1] - ray->org[1]) / ray->dir[1];
    }
  if (ray->dir[2] < 0.0)
    {
      inc[2] = -1;
      t_max[2] = (i[2] * mesh->voxel_size - mesh->half_dim[2] - ray->org[2]) / ray->dir[2];
    }
  else
    {
      inc[2] = 1;
      t_max[2] = ((i[2] + 1) * mesh->voxel_size - mesh->half_dim[2] - ray->org[2]) / ray->dir[2];
    }
  for (l = 0; l < 3; l++)
    t_delta[l] = inc[l] * mesh->voxel_size / ray->dir[l];
  
  for (;;)
    {
      if (ray->t_near > ray->t_max) return !SUCCESS;
      
      ray->t_far = fmin(t_max[0], fmin(t_max[1], t_max[2]));
      
      voxel = &mesh->voxel[ (i[0]*mesh->voxels[1]+i[1])*mesh->voxels[2]+i[2] ];
      
      if (voxel->triangles > 0 &&
	  rtVoxelVsRay (voxel, ray, hit, mesh))
	{
	  return SUCCESS;
	}
      if (t_max[0] < t_max[1])
	{
	  if (t_max[0] < t_max[2])
	    {
	      i[0] += inc[0];
	      
	      if (i[0] < 0 || i[0] >= mesh->voxels[0]) return !SUCCESS;
	      
	      t_max[0] += t_delta[0];
	    }
	  else
	    {
	      i[2] += inc[2];
	      
	      if (i[2] < 0 || i[2] >= mesh->voxels[2]) return !SUCCESS;
	      
	      t_max[2] += t_delta[2];
	    }
	}
      else
	{
	  if (t_max[1] < t_max[2])
	    {
	      i[1] += inc[1];
	      
	      if (i[1] < 0 || i[1] >= mesh->voxels[1]) return !SUCCESS;
	      
	      t_max[1] += t_delta[1];
	    }
	  else
	    {
	      i[2] += inc[2];
	      
	      if (i[2] < 0 || i[2] >= mesh->voxels[2]) return !SUCCESS;
	      
	      t_max[2] += t_delta[2];
	    }
	}
      ray->t_near = ray->t_far;
    }
}


int rtTracePrimaryRay (Hit *first_hit, Vis *vis)
{
  double screen_dir_x[3], screen_dir_y[3];
  double screen_vtx[3];
  double aabb[6];
  
  int l;
  
  Ray ray;
  
  
  Rotate (screen.dim[0], 0.0, 0.0,
	  viewpoint.sin_longitude, viewpoint.cos_longitude,
	  viewpoint.sin_latitude, viewpoint.cos_latitude,
	  &screen_dir_x[0], &screen_dir_x[1], &screen_dir_x[2]);
  
  Rotate (0.0, screen.dim[1], 0.0,
	  viewpoint.sin_longitude, viewpoint.cos_longitude,
	  viewpoint.sin_latitude, viewpoint.cos_latitude,
	  &screen_dir_y[0], &screen_dir_y[1], &screen_dir_y[2]);
  
  for (l = 0; l < 3; l++)
    {
      screen_vtx[l] = (-viewpoint.dist / vis->viewpoint_radius) * viewpoint.pos[l] -
	screen_dir_x[l] - screen_dir_y[l];
      
      screen_dir_x[l] *= (2.0 / screen.pixels[0]);
      screen_dir_y[l] *= (2.0 / screen.pixels[1]);
    }
  for (l = 0; l < 3; l++)
    ray.org[l] = viewpoint.pos[l];
  
  for (l = 0; l < 3; l++)
    ray.dir[l] = screen_vtx[l] + vis->mouse.x[0]*screen_dir_x[l] + vis->mouse.x[1]*screen_dir_y[l];
  
  ray.t_max = 1.0e+30;
  
  aabb[0] =  vis->mesh.half_dim[0];
  aabb[1] = -vis->mesh.half_dim[0];
  aabb[2] =  vis->mesh.half_dim[1];
  aabb[3] = -vis->mesh.half_dim[1];
  aabb[4] =  vis->mesh.half_dim[2];
  aabb[5] = -vis->mesh.half_dim[2];
  
  rtAABBvsRay (aabb, &ray);
  
  first_hit->previous_triangle_id = -1;
  
  if (ray.t_near > ray.t_far ||
      rtTraceRay (&ray, first_hit, &vis->mesh) != SUCCESS)
    {
      return !SUCCESS;
    }
  return SUCCESS;
}


int rtTraceSecondaryRay (Hit *first_hit, Hit *second_hit, Vis *vis)
{
  double primary_dir[3];
  
  int l;
  
  Ray ray;
  
  
  for (l = 0; l < 3; l++)
    ray.org[l] = first_hit->pos[l];
  
  for (l = 0; l < 3; l++)
    primary_dir[l] = first_hit->pos[l] - viewpoint.pos[l];
  
  for (l = 0; l < 3; l++)
    ray.dir[l] = vis->mesh.triangle[ first_hit->triangle_id ].nor[l];
  
  if (ScalarProd (primary_dir, ray.dir) < 0.0)
    {
      for (l = 0; l < 3; l++)
	ray.dir[l] = -ray.dir[l];
    }
  ray.t_max = 1.0e+30;
  ray.t_near = 0.0;
  
  second_hit->previous_triangle_id = first_hit->triangle_id;
  
  if (rtTraceRay (&ray, second_hit, &vis->mesh) != SUCCESS)
    {
      return !SUCCESS;
    }
  return SUCCESS;
}


void rtEndRayTracing (Mesh *mesh)
{
  for (int n = 0; n < mesh->voxels[3]; n++)
    {
      if (mesh->voxel[n].triangle_id != NULL)
	{
	  free(mesh->voxel[n].triangle_id);
	  mesh->voxel[n].triangle_id = NULL;
	}
    }
  free(mesh->voxel);
}


#endif // MESH
