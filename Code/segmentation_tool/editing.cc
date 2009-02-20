#include "editing.h"


void editTriangleCenter (Triangle *t_p, float c[3])
{
  for (int l = 0; l < 3; l++)
    c[l] = (1.0F/3.0F) * (t_p->v[0].pos[l]+t_p->v[1].pos[l]+t_p->v[2].pos[l]);
}


void editTriangleNormal (Triangle *t_p, float nor[3])
{
  float dx1[3], dx2[3];
  float temp;
  
  
  for (int l = 0; l < 3; l++)
    {
      dx1[l] = t_p->v[1].pos[l] - t_p->v[0].pos[l];
      dx2[l] = t_p->v[2].pos[l] - t_p->v[0].pos[l];
    }
  nor[0] = dx2[1] * dx1[2] - dx2[2] * dx1[1];
  nor[1] = dx2[2] * dx1[0] - dx2[0] * dx1[2];
  nor[2] = dx2[0] * dx1[1] - dx2[1] * dx1[0];
  
  temp = 1.F / fmaxf(1.0e-30F, sqrtf(nor[0]*nor[0] + nor[1]*nor[1] + nor[2]*nor[2]));
  
  nor[0] *= temp;
  nor[1] *= temp;
  nor[2] *= temp;
}


float editTriangleRadius (Triangle *t_p)
{
  float dx[3][3];
  float mx[3][3];
  float c[3], r2[3];
  
  
  for (int l = 0; l < 3; l++)
    {
      mx[0][l] = 0.5F * (t_p->v[1].pos[l] + t_p->v[0].pos[l]);
      mx[1][l] = 0.5F * (t_p->v[2].pos[l] + t_p->v[0].pos[l]);
      mx[2][l] = 0.5F * (t_p->v[2].pos[l] + t_p->v[1].pos[l]);
    }
  editTriangleCenter (t_p, c);
  
  for (int m = 0; m < 3; m++)
    for (int l = 0; l < 3; l++)
      dx[m][l] = mx[m][l] - c[l];
  
  for (int m = 0; m < 3; m++)
    r2[m] = dx[m][0]*dx[m][0] + dx[m][1]*dx[m][1] + dx[m][2]*dx[m][2];
  
  if (r2[0] < r2[1] && r2[0] < r2[2])
    {
      return sqrtf(r2[0]);
    }
  else if (r2[1] < r2[0] && r2[1] < r2[2])
    {
      return sqrtf(r2[1]);
    }
  else
    {
      return sqrtf(r2[2]);
    }
}


void editCalculateTriangleData (Triangle *t_p)
{
  float longitude, latitude;
  float r;
  
  
  editTriangleCenter (t_p, t_p->pos);
  editTriangleNormal (t_p, t_p->nor);
  
  longitude = atan2f(t_p->nor[0], t_p->nor[2]);
  latitude  = atan2f(t_p->nor[1], sqrtf(t_p->nor[0]*t_p->nor[0] + t_p->nor[2]*t_p->nor[2]));
  
  t_p->d.sin_longitude = sinf(longitude);
  t_p->d.cos_longitude = cosf(longitude);
  t_p->d.sin_latitude  = sinf(latitude);
  t_p->d.cos_latitude  = cosf(latitude);
  
  r = editTriangleRadius (t_p);
  
  t_p->d.r2 = r * r;
}


void editDeleteTriangle (int b_id, int t_id, Vis *vis)
{
  Boundary *boundary_p = &vis->boundary[ b_id ];
  
  if (t_id == --boundary_p->triangles) return;
  
  memcpy (&boundary_p->triangle[ t_id ],
	  &boundary_p->triangle[ boundary_p->triangles ], sizeof(Triangle));
}


void editInvertTriangleNormal (int b_id, int t_id, Vis *vis)
{
  Triangle *t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  t_p->normal_sign = -t_p->normal_sign;
}


void editChangeTrianglePars (int b_id, int t_id, float dp_avg, float dp_amp, float dp_phs, Vis *vis)
{
  Triangle *t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  t_p->pressure_avg += dp_avg;
  t_p->pressure_avg = fmaxf(0.0, 0.1F * (int)(10.F * t_p->pressure_avg));
  
  t_p->pressure_amp += dp_amp;
  t_p->pressure_amp = fminf(t_p->pressure_avg, fmaxf(0.0F, 0.1F * (int)(10.F * t_p->pressure_amp)));
  
  t_p->pressure_phs += dp_phs;
  t_p->pressure_phs = fmaxf(0.0, 0.1F * (int)(10.F * t_p->pressure_phs));
}


void editRescaleSystem (Vis *vis)
{
  int l;
  
  
  vis->output_voxels[0] = vis->input_voxels[0] * vis->res_factor;
  vis->output_voxels[1] = vis->input_voxels[1] * vis->res_factor;
  vis->output_voxels[2] = (int)(vis->input_voxels[2] * (vis->slice_size / vis->pixel_size) *
				(float)vis->output_voxels[0]/(float)vis->input_voxels[0]);
  
  vis->scale[0] = vis->res_factor;
  vis->scale[1] = vis->res_factor;
  vis->scale[2] = vis->res_factor * vis->slice_size / vis->pixel_size;
  
  for (l = 0; l < 3; l++)
    vis->inv_scale[l] = 1.0 / vis->scale[l];
  
  for (l = 0; l < 3; l++)
    vis->blocks[l] = vis->output_voxels[l] >> SHIFT;
  
  for (l = 0; l < 3; l++)
    if ((vis->blocks[l] << SHIFT) < vis->output_voxels[l]) ++vis->blocks[l];
  
  for (l = 0; l < 3; l++)
    vis->sites[l] = vis->blocks[l] * BLOCK_SIZE;
  
  for (l = 0; l < 3; l++)
    vis->dim[l] = vis->output_voxels[l];
  
  for (l = 0; l < 3; l++)
    vis->half_dim[l] = 0.5F * vis->dim[l];
  
  vis->system_size = fmaxf(vis->dim[0], fmaxf(vis->dim[1], vis->dim[2]));
  
  int temp = vis->tot_blocks;
  
  vis->tot_blocks = vis->blocks[0] * vis->blocks[1] * vis->blocks[2];
  
  if (temp < vis->tot_blocks)
    {
      vis->block = (Block *)realloc(vis->block, sizeof(Block) * vis->tot_blocks);
    }
}


void editRescaleTriangles (float scale, Vis *vis)
{
  for (int n = 0; n < BOUNDARIES; n++)
    for (int m = 0; m < vis->boundary[ n ].triangles; m++)
      {
	for (int l = 0; l < 3; l++)
	  for (int k = 0; k < 3; k++)
	    vis->boundary[n].triangle[m].v[l].pos[k] *= scale;
	
	editCalculateTriangleData (&vis->boundary[n].triangle[m]);
      }
}


void editProjectBoundariesToScreenVoxels (Vis *vis)
{
  float x[3];
  
  int voxel[2];
  
  ScreenVoxel *voxel_p;
  
  
  for (int n = 0; n < vis->screen_voxels * vis->screen_voxels; n++)
    {
      vis->screen_voxel[ n ].z[VERTEX] = 0.F;
    }
  for (int n = 0; n < BOUNDARIES; n++)
    {
      for (int m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  Triangle *t_p = &vis->boundary[ n ].triangle[ m ];
	  
	  for (int l = 0; l < 3; l++)
	    {
	      visProject (t_p->v[l].pos, x);
	      
	      voxel[0] = (int)((vis->screen_voxels/(2.0F*screen.dim[0])) * (x[0]+screen.dim[0]));
	      voxel[1] = (int)((vis->screen_voxels/(2.0F*screen.dim[1])) * (x[1]+screen.dim[1]));
	      
	      if (voxel[0] < 0 || voxel[0] >= vis->screen_voxels ||
		  voxel[1] < 0 || voxel[1] >= vis->screen_voxels)
		{
		  continue;
		}
	      voxel_p = &vis->screen_voxel[ voxel[0]*vis->screen_voxels + voxel[1] ];
	      
	      if (x[2] > voxel_p->z[VERTEX])
		{
		  voxel_p->b_id = n;
		  voxel_p->t_id = m;
		  voxel_p->v_id = l;
		  voxel_p->z[VERTEX] = x[2];
		}
	    }
	}
    }
}


void editMoveTriangleVertexWithMouse (int x0[3], Vis *vis)
{
  float x1[3];
  float x2[3];
  
  int i, j, k;
  
  
  x1[0] = screen.dim[0] * (-1.0 + (x0[0]<<1) / (float)vis->viewport_pixels[0]);
  x1[1] = screen.dim[1] * (-1.0 + (x0[1]<<1) / (float)vis->viewport_pixels[1]);
  x1[2] = x0[2];
  
  visAntiProject (x1, x2);
  
  i = vis->mouse.b_id;
  j = vis->mouse.t_id;
  k = vis->mouse.v_id;
  
  for (int l = 0; l < 3; l++)
    vis->boundary[i].triangle[j].v[k].pos[l] = x2[l];
  
  editCalculateTriangleData (&vis->boundary[i].triangle[j]);
}


void editRotateTriangleWithMouse (float dx0, float dy0, int b_id, int t_id, Vis *vis)
{
  float c[3], nx[3], x[3];
  float dx1[3], dx2[3], dx3[3];
  float longitude[2], latitude[2];
  
  int i, l;
  
  Triangle *t_p;
  
  
  t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  editTriangleCenter (t_p, c);
  editTriangleNormal (t_p, nx);
  
  longitude[0] = atan2f(nx[0], nx[2]);
  latitude[0]  = atan2f(nx[1], sqrtf(nx[0]*nx[0] + nx[2]*nx[2]));
  
  x[0] = dx0;
  x[1] = dy0;
  x[2] = 1.0F;
  
  Rotate (x, longitude[0], latitude[0], dx1);
  
  longitude[1] = atan2f(dx1[0], dx1[2]);
  latitude[1]  = atan2f(dx1[1], sqrtf(dx1[0]*dx1[0] + dx1[2]*dx1[2]));
  
  for (i = 0; i < 3; i++)
    {
      for (l = 0; l < 3; l++)
	dx1[l] = t_p->v[i].pos[l] - c[l];
      
      AntiRotate (dx1, longitude[0], latitude[0], dx2);
      Rotate (dx2, longitude[1], latitude[1], dx3);
      
      for (l = 0; l < 3; l++)
	t_p->v[i].pos[l] = dx3[l] + c[l];
    }
  editCalculateTriangleData (t_p);
}


void editScaleTriangleWithMouse (float scaling_factor, int b_id, int t_id, Vis *vis)
{
  float c[3];
  
  int i, l;
  
  Triangle *t_p;
  
  
  t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  editTriangleCenter (t_p, c);
  
  for (i = 0; i < 3; i++)
    for (l = 0; l < 3; l++)
      t_p->v[i].pos[l] = (1.0+scaling_factor) * (t_p->v[i].pos[l]-c[l]) + c[l];
  
  editCalculateTriangleData (t_p);
}

/*
void editRotateViewpointWithMouse (float dx0, float dy0, Vis *vis)
{
  float dx1[3], x[3];
  
  
  x[0] = dx0;
  x[1] = dy0;
  x[2] = 1.0F;
  
  Rotate (x, vis->longitude * DEG_TO_RAD, vis->latitude * DEG_TO_RAD, dx1);
  
  vis->longitude = atan2f(dx1[0], dx1[2]) / DEG_TO_RAD;
  vis->latitude  = atan2f(dx1[1], sqrtf(dx1[0]*dx1[0] + dx1[2]*dx1[2])) / DEG_TO_RAD;
  
  visProjection (vis);
}
*/

void editRotateViewpointWithMouse (float dx, float dy, Vis *vis)
{
  vis->longitude += dx / DEG_TO_RAD;
  vis->latitude  += dy / DEG_TO_RAD;
  
  visProjection (vis);
}


void editMouseFunction (int button, int state, int x, int y, Vis *vis)
{
  int voxel[2];
  
  
  vis->mouse.state = !ACTIVE;
  
  if (button != GLUT_LEFT_BUTTON) return;
  
  y = vis->viewport_pixels[1] - y - 1;
  
  voxel[0] = (x * vis->screen_voxels) / vis->viewport_pixels[0];
  voxel[1] = (y * vis->screen_voxels) / vis->viewport_pixels[1];
  
  if (voxel[0] < 0 || voxel[0] >= vis->screen_voxels ||
      voxel[1] < 0 || voxel[1] >= vis->screen_voxels)
    {
      vis->mouse.b_id = -1;
      return;
    }
  if (state == GLUT_DOWN)
    {
      vis->mouse.state = ACTIVE;
    }
  else if (state == GLUT_UP)
    {
      vis->mouse.state = !ACTIVE;
    }
  vis->mouse.x[0] = x;
  vis->mouse.x[1] = y;
}


void editMotionFunction (int x, int y, Vis *vis)
{
  if (vis->mouse.state == !ACTIVE) return;
  
  y = vis->viewport_pixels[1] - y - 1;
  
  int mouse_dx = x - vis->mouse.x[0];
  int mouse_dy = y - vis->mouse.x[1];
  
  if (vis->menu.option == CHANGE_THRESHOLD)
    {
      vis->mouse.dy = mouse_dy;
    }
  else if (vis->mode == 0)
    {
      if (vis->menu.option == CHANGE_SLICE)
	{
	  vis->selected_voxel[2] += mouse_dy;
	  vis->selected_voxel[2] = max(0, min(vis->input_voxels[2]-1, vis->selected_voxel[2]));
	}
    }
  else
    {
      if (vis->menu.option == ZOOM_SCENE)
	{
	  vis->zoom *= 1.0 + (float)mouse_dy / vis->viewport_pixels[1];
	  visProjection (vis);
	}
      else if (vis->menu.option == ROTATE_SCENE)
	{
	  editRotateViewpointWithMouse (-(float)mouse_dx / vis->viewport_pixels[0],
					-(float)mouse_dy / vis->viewport_pixels[1], vis);
	}
      else if (vis->menu.option == ZOOM_BOUNDARY && vis->mouse.b_id >= 0)
	{
	  editScaleTriangleWithMouse ((float)mouse_dy / vis->viewport_pixels[1],
				      vis->mouse.b_id,
				      vis->mouse.t_id,
				      vis);
	}
      else if (vis->menu.option == ROTATE_BOUNDARY && vis->mouse.b_id >= 0)
	{
	  editRotateTriangleWithMouse ((float)mouse_dx / vis->viewport_pixels[0],
				       (float)mouse_dy / vis->viewport_pixels[1],
				       vis->mouse.b_id,
				       vis->mouse.t_id,
				       vis);
	}
      else if ((vis->menu.option & (CHANGE_PRESSURE_AMPLITUDE|CHANGE_MEAN_PRESSURE|CHANGE_PRESSURE_PHASE)) &&
	       vis->mouse.b_id >= 0 &&
	       (vis->mouse.b_id == INLET_BOUNDARY || vis->mouse.b_id == OUTLET_BOUNDARY))
	{
	  float dy = 10.F * (float)mouse_dy / vis->viewport_pixels[1];
	  
	  if (vis->menu.option & CHANGE_MEAN_PRESSURE)
	    {
	      editChangeTrianglePars (vis->mouse.b_id, vis->mouse.t_id, dy, 0.F, 0.F, vis);
	    }
	  else if (vis->menu.option & CHANGE_PRESSURE_AMPLITUDE)
	    {
	      editChangeTrianglePars (vis->mouse.b_id, vis->mouse.t_id, 0.F, dy, 0.F, vis);
	    }
	  else
	    {
	      editChangeTrianglePars (vis->mouse.b_id, vis->mouse.t_id, 0.F, 0.F, dy, vis);
	    }
	}
    }
  if (vis->mode == 2)
    {
      vis->mode = 1;
    }
  vis->mouse.x[0] = x;
  vis->mouse.x[1] = y;
}


void editPassiveMotionFunction (int x, int y, Vis *vis)
{
  int voxel[2];
  
  ScreenVoxel *voxel_p;
  
  
  y = vis->viewport_pixels[1] - y - 1;
  
  vis->mouse.x[0] = x;
  vis->mouse.x[1] = y;
  
  if (vis->mode == 0) return;
  
  voxel[0] = (x * vis->screen_voxels) / vis->viewport_pixels[0];
  voxel[1] = (y * vis->screen_voxels) / vis->viewport_pixels[1];
  
  if (voxel[0] < 0 || voxel[0] >= vis->screen_voxels ||
      voxel[1] < 0 || voxel[1] >= vis->screen_voxels)
    {
      vis->mouse.b_id = -1;
      return;
    }
  if (vis->mouse.b_id >= 0) return;
  
  editProjectBoundariesToScreenVoxels (vis);
  
  voxel_p = visScreenVoxelPtr (vis->mouse.x, vis);
  
  if (voxel_p->z[VERTEX] > EPSILON)
    {
      vis->mouse.t_id = voxel_p->t_id;
      vis->mouse.b_id = voxel_p->b_id;
      vis->mouse.v_id = voxel_p->v_id;
    }
}


void GLUTCALLBACK MouseFunction (int button, int state, int x, int y)
{
  editMouseFunction (button, state, x, y, &vis);
}


void GLUTCALLBACK MotionFunction (int x, int y)
{
  editMotionFunction (x, y, &vis);
}


void GLUTCALLBACK PassiveMotionFunction (int x, int y)
{
  editPassiveMotionFunction (x, y, &vis);
}


void GLUTCALLBACK Reshape (GLsizei w, GLsizei h)
{
  // the window is reshaped if necessary
  
  vis.ortho[0] *= (float)w / (float)vis.viewport_pixels[0];
  vis.ortho[1] *= (float)h / (float)vis.viewport_pixels[1];
  
  vis.viewport_pixels[0] = w;
  vis.viewport_pixels[1] = h;
  
  glViewport(0, 0, w, h);
  
  visProjection (&vis);
}
