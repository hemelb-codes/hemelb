#include "vis.h"


void visColorPalette (int iters, int res_factor, double col[3])
{
  double t;
  
  
  iters = (iters / res_factor) % 200;
  
  if (iters > 100) iters = 200 - iters;
  
  t = 0.01 * (double)iters;
  
  if (t > 1.0)
    {
      col[0] = 1.0;
      col[1] = 0.0;
      col[2] = 0.0;
    }
  else if (t > 0.0 && t <= 0.25)
    {
      col[0] = 0.0;
      col[1] = 4.0 * t;
      col[2] = 1.0;
    }
  else if (t > 0.25 && t <= 0.5)
    {
      col[0] = 0.0;
      col[1] = 1.0;
      col[2] = 2.0 - 4.0 * t;
    }
  else if (t > 0.5 && t <= 0.75)
    {
      col[0] = 4.0 * (t - 0.5);
      col[1] = 1.0;
      col[2] = 0.0;
    }
  else if (t > 0.75 && t <= 1.0)
    {
      col[0] = 1.0;
      col[1] = 4.0 - 4.0 * t;
      col[2] = 0.0;
    }
  else
    {
      col[0] = 0.0;
      col[1] = 0.0;
      col[2] = 1.0;
    }
}


ScreenVoxel *visScreenVoxelPtr (short int x[2], Vis *vis)
{
  int voxel[2];
  
  for (int l = 0; l < 2; l++)
    voxel[l] = (x[l] * vis->screen_voxels) / vis->viewport_pixels[l];
  
  return &vis->screen_voxel[ voxel[0]*vis->screen_voxels+voxel[1] ];
}


void visOpenWindow (int pixels_x, int pixels_y)
{
  screen.pixels[0] = pixels_x;
  screen.pixels[1] = pixels_y;
  
  glutInitDisplayMode (GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutInitWindowPosition (0, 0);
  glutInitWindowSize (pixels_x, pixels_y);
  
  glutCreateWindow (" ");
  
  glEnable (GL_DEPTH_TEST);
  glDisable (GL_BLEND);
  glShadeModel (GL_SMOOTH);
  glDisable(GL_DITHER);
  glDisable (GL_LIGHTING);
  
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}


void visProjection (Vis *vis)
{
  double temp;
  
  
  for (int l = 0; l < 3; l++)
    screen.col[l] = vis->background[l];
  
  screen.dim[0] = vis->ortho[0] / vis->zoom;
  screen.dim[1] = vis->ortho[1] / vis->zoom;
  
  screen.pixels[0] = vis->viewport_pixels[0];
  screen.pixels[1] = vis->viewport_pixels[1];
  
  viewpoint.sin_longitude = sin(vis->longitude * DEG_TO_RAD);
  viewpoint.cos_longitude = cos(vis->longitude * DEG_TO_RAD);
  
  viewpoint.sin_latitude = sin(vis->latitude * DEG_TO_RAD);
  viewpoint.cos_latitude = cos(vis->latitude * DEG_TO_RAD);
  
  temp = vis->viewpoint_radius * viewpoint.cos_latitude;
  
  viewpoint.pos[0] = temp * viewpoint.sin_longitude;
  viewpoint.pos[1] = vis->viewpoint_radius * viewpoint.sin_latitude;
  viewpoint.pos[2] = temp * viewpoint.cos_longitude;
  
  for (int l = 0; l < 3; l++)
    viewpoint.pos[l] += vis->scene_center[l];
  
  viewpoint.dist = vis->viewport_radius / 2.0;
  
  temp = vis->viewport_radius / vis->viewpoint_radius;
  
  for (int l = 0; l < 3; l++)
    screen.ctr[l] = viewpoint.pos[l] + temp*(vis->scene_center[l] - viewpoint.pos[l]);
  
  screen.zoom = vis->zoom;
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  
  glFrustum (-screen.dim[0], screen.dim[0],
	     -screen.dim[1], screen.dim[1],
	     viewpoint.dist, 1.0e+30);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  if ((int)((fabs(vis->latitude) + 90.0) / 180.0) % 2 == 0)
    {
      gluLookAt (viewpoint.pos[0], viewpoint.pos[1], viewpoint.pos[2], 
		 screen.ctr[0], screen.ctr[1], screen.ctr[2],
		 0.0, 1.0, 0.0);
    }
  else
    {
      gluLookAt (viewpoint.pos[0], viewpoint.pos[1], viewpoint.pos[2], 
		 screen.ctr[0], screen.ctr[1], screen.ctr[2],
		 0.0, -1.0, 0.0);
    }
  glClearColor (screen.col[0], screen.col[1], screen.col[2], 0.0);
}


void visProject (double px1[3], double px2[3])
{
  double x[3];
  
  
  for (int l = 0; l < 3; l++)
    x[l] = px1[l] - viewpoint.pos[l];
  
  AntiRotate (x[0], x[1], x[2],
	      viewpoint.sin_longitude, viewpoint.cos_longitude,
	      viewpoint.sin_latitude, viewpoint.cos_latitude,
	      &px2[0], &px2[1], &px2[2]);
  
  double temp = viewpoint.dist / (px2[2] = -px2[2]);
  
  px2[0] *= temp;
  px2[1] *= temp;
}


void visAntiProject (double px1[3], double px2[3])
{
  double x[3];
  double temp;
  
  
  temp = px1[2] / viewpoint.dist;
  
  x[0] = temp * px1[0];
  x[1] = temp * px1[1];
  x[2] = -px1[2];
  
  Rotate (x[0], x[1], x[2],
	  viewpoint.sin_longitude, viewpoint.cos_longitude,
	  viewpoint.cos_latitude, viewpoint.cos_latitude,
	  &px2[0], &px2[1], &px2[2]);
  
  for (int l = 0; l < 3; l++)
    px2[l] += viewpoint.pos[l];
}


void visCalculateSceneCenter (Vis *vis)
{
  int block_min[3], block_max[3];
  int b[3];
  
  double scale;
  
  
  for (int l = 0; l < 3; l++)
    {
      block_min[l] = +100000000;
      block_max[l] = -100000000;
    }
  int n = -1;
  
  for (b[0] = 0; b[0] < vis->blocks[0]; b[0]++)
    for (b[1] = 0; b[1] < vis->blocks[1]; b[1]++)
      for (b[2] = 0; b[2] < vis->blocks[2]; b[2]++)
	{
	  if (vis->block[++n].site == NULL) continue;
	  
	  for (int l = 0; l < 3; l++)
	    {
	      block_min[l] = min(block_min[l], b[l]);
	      block_max[l] = max(block_max[l], b[l]);
	    }
	}
#ifndef MESH
  scale = 1.0 / vis->res_factor;
#else
  scale = vis->mesh.voxel_size / vis->res_factor;
#endif
  for (int l = 0; l < 3; l++)
    vis->scene_center[l] = (BLOCK_SIZE * 0.5 * (block_min[l] + block_max[l])) * scale - vis->half_dim[l];
}



void visCreateCubeDisplayList (Vis *vis)
{
#ifndef MESH
  double scale = 0.5 / vis->res_factor;
#else
  double scale = 0.5 * vis->mesh.voxel_size / vis->res_factor;
#endif
  glNewList (1, GL_COMPILE);
  
  glBegin (GL_QUAD_STRIP);
  glVertex3d (-scale, -scale, -scale);
  glVertex3d (-scale, +scale, -scale);
  glVertex3d (-scale, -scale, +scale);
  glVertex3d (-scale, +scale, +scale);
  glVertex3d (+scale, -scale, +scale);
  glVertex3d (+scale, +scale, +scale);
  glVertex3d (+scale, -scale, -scale);
  glVertex3d (+scale, +scale, -scale);
  glVertex3d (-scale, -scale, -scale);
  glVertex3d (-scale, +scale, -scale);
  glEnd ();
  
  glBegin (GL_QUADS);
  glVertex3d (-scale, -scale, -scale);
  glVertex3d (-scale, -scale, +scale);
  glVertex3d (+scale, -scale, +scale);
  glVertex3d (+scale, -scale, -scale);
  glVertex3d (-scale, +scale, -scale);
  glVertex3d (-scale, +scale, +scale);
  glVertex3d (+scale, +scale, +scale);
  glVertex3d (+scale, +scale, -scale);
  glEnd ();
  
  glEndList ();
}


void visDeleteCubeDisplayList (void)
{
  glDeleteLists (1, 1);
}


void visVisualiseString (double r, double g, double b, int x, int y, char *string, void *font)
{
  glColor3d (r, g, b);
  glWindowPos2i (x, y);
  
  for (int i = 0; i < (int)strlen(string); i++)
    {
      glutBitmapCharacter (font, string[i]);
    }
  glEnd ();
}


void visVisualiseTrianglePars (int b_id, int t_id, Vis *vis)
{
  char pars_string[256];
  
  BoundaryTriangle *t_p;
  
  
  t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  sprintf (pars_string, "pressure = %.1f + %.1f cos(w t + phase)  mmHg, phase = %.1f deg",
	   t_p->pressure_avg, t_p->pressure_amp, t_p->pressure_phs);
  
  visVisualiseString (0.0, 0.0, 0.0, 5, 5,
		      pars_string, GLUT_BITMAP_HELVETICA_12);
}


#ifndef MESH
void visVisualiseSiteData (short int site[3], Vis *vis)
#else
void visVisualiseSiteData (Vis *vis)
#endif
{
  double diameter;
  
  char pars_string[256];
  
  
#ifndef MESH
  double nor[3];
  
  if (segEstimateNormal (site, nor, vis) != SUCCESS)
    {
      return;
    }
  segEstimateDiameter (site, nor, &diameter, vis);
  
  if (diameter >= 0.0)
    {
      diameter *= vis->pixel_size / vis->res_factor;
      sprintf (pars_string, "Lattice coords = (%i, %i, %i)  Diameter = %.3f mm",
	       site[0], site[1], site[2], diameter);
    }
  else
    {
      sprintf (pars_string, "Lattice coords = (%i, %i, %i)  Diameter = NaN",
	       site[0], site[1], site[2]);
    }
#else // MESH
  Hit first_hit, second_hit;
  
  if (segEstimateExtrema (&first_hit, &second_hit, vis) == SUCCESS)
    {
      segEstimateDiameter (&diameter, vis);
      sprintf (pars_string, "Diameter = %.3f mm", diameter);
      
      visVisualiseHitData (&first_hit, &second_hit, vis);
    }
  else
    {
      sprintf (pars_string, "Diameter = NaN");
    }
#endif // MESH
  visVisualiseString (0.0, 0.0, 0.0, 5, 5,
		      pars_string, GLUT_BITMAP_HELVETICA_12);
}


void visVisualiseVisData (Vis *vis)
{
  char data_string[256];
  
  
  if (vis->mode == 0)
    {
#ifndef MESH
      sprintf (data_string, "Pixel size (mm) = %.3f  Slice thickness (mm) = %.3f  Threshold = %.1f  ",
	       vis->pixel_size, vis->slice_size, vis->selected_grey);
#else
      sprintf (data_string, "Mesh voxel size (mm) = %e", vis->mesh.voxel_size);
#endif
      visVisualiseString (0.0, 1.0, 0.0,
			  5, vis->viewport_pixels[1]-20,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "Resolution enhancement = %i  #sites = %i  #sup. sites = %i  Seg. time = %.3f s",
	       vis->res_factor, vis->tot_sites, vis->coords[C], vis->segmentation_time);
      
      visVisualiseString (0.0, 1.0, 0.0, 5, vis->viewport_pixels[1]-40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
    }
  else
    {
      sprintf (data_string, "longitude = %.1f", vis->longitude);
      
      visVisualiseString (0.0, 0.0, 0.0, vis->viewport_pixels[0]-110, 55,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "latitude = %.1f", vis->latitude);
      
      visVisualiseString (0.0, 0.0, 0.0, vis->viewport_pixels[0]-110, 40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "zoom = %.1f", vis->zoom);
      
      visVisualiseString (0.0, 0.0, 0.0, vis->viewport_pixels[0]-110, 25,
			  data_string, GLUT_BITMAP_HELVETICA_12);
#ifndef MESH
      sprintf (data_string, "Pixel size (mm) = %.3f  Slice thickness (mm) = %.3f  Threshold = %.1f  ",
	       vis->pixel_size, vis->slice_size, vis->selected_grey);
#else
      sprintf (data_string, "Mesh voxel size (mm) = %.3f",
	       vis->mesh.voxel_size);
#endif
      visVisualiseString (0.0, 0.0, 0.0,
			  5, vis->viewport_pixels[1]-20,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "Resolution enhancement = %i  #sites = %i  #sup. sites = %i  Seg. time = %.3f s",
	       vis->res_factor, vis->tot_sites, vis->coords[C], vis->segmentation_time);
      
      visVisualiseString (0.0, 0.0, 0.0, 5, vis->viewport_pixels[1]-40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "FPS = %.2f", vis->fps);
      
      visVisualiseString (0.0, 0.0, 0.0,
			  vis->viewport_pixels[0]-110, 5,
			  data_string, GLUT_BITMAP_HELVETICA_12);
    }
}


void visVisualiseTriangles (Vis *vis)
{
  double x1[3], x2[3];
  double length;
  
  BoundaryTriangle *t_p;
  
  
  glLineWidth (2.0);
  
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      if (n == INLET_BOUNDARY)
	{
	  glColor3d (0.5, 0.5, 0.5);
	}
      else if (n == OUTLET_BOUNDARY)
	{
	  glColor3d (0.0, 0.5, 0.0);
	}
      else if (n == WALL_BOUNDARY)
	{
	  glColor3d (0.0, 0.0, 0.5);
	}
      for (unsigned int m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  if (m == vis->mouse.t_id && n == vis->mouse.b_id)
	    {
	      glBegin (GL_TRIANGLES);
	    }
	  else
	    {
	      glBegin (GL_LINE_LOOP);
	    }
	  t_p = &vis->boundary[ n ].triangle[ m ];
	  
	  for (int i = 0; i < 3; i++)
	    glVertex3dv (t_p->v[i].pos);
	  
	  glEnd ();
	  
	  if (n != INLET_BOUNDARY) continue;
	  
	  glBegin (GL_LINES);
	  
	  length = sqrt(t_p->d.r2);
	  
	  for (int l = 0; l < 3; l++)
	    x2[l] = (x1[l] = t_p->pos[l]) + t_p->normal_sign * t_p->nor[l] * length;
	  
	  glVertex3dv (x1);
	  glVertex3dv (x2);
	  
	  glEnd ();
	}
    }
}


void visVisualiseDiscs (Vis *vis)
{
  double v[2*37];
  double x1[3], x2[3];
  
  BoundaryTriangle *t_p;
  
  
  for (int i = 0; i <= 35; i++)
    {
      v[ 2*i   ] = sin(i * 10 * DEG_TO_RAD);
      v[ 2*i+1 ] = cos(i * 10 * DEG_TO_RAD);
    }
  v[ 2*36   ] = v[ 0 ];
  v[ 2*36+1 ] = v[ 1 ];
  
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      if (n == INLET_BOUNDARY)
	{
	  glColor3d (0.5, 0.5, 0.5);
	}
      else if (n == OUTLET_BOUNDARY)
	{
	  glColor3d (0.0, 0.5, 0.0);
	}
      else if (n == WALL_BOUNDARY)
	{
	  glColor3d (0.0, 0.0, 0.5);
	}
      for (unsigned int m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  glBegin (GL_TRIANGLE_FAN);
	  
	  t_p = &vis->boundary[ n ].triangle[ m ];
	  
	  glVertex3dv (t_p->pos);
	  
	  for (int i = 0; i <= 36; i++)
	    {
	      x1[0] = sqrt(t_p->d.r2) * v[ 2*i   ];
	      x1[1] = sqrt(t_p->d.r2) * v[ 2*i+1 ];
	      x1[2] = 0.0;
	      
	      Rotate (x1[0], x1[1], x1[2],
		      t_p->d.sin_longitude, t_p->d.cos_longitude,
		      t_p->d.sin_latitude,  t_p->d.cos_latitude,
		      &x2[0], &x2[1], &x2[2]);
	      
	      for (int l = 0; l < 3; l++)
		x1[l] = x2[l] + t_p->pos[l];
	      
	      glVertex3dv (x1);
	    }
	  glEnd ();
	}
    }
}


void visVisualiseActiveBoundaryVoxel (Vis *vis)
{
  double x1[2], x2[2];
  
  int matrix_mode;
  
  
  glGetIntegerv (GL_MATRIX_MODE, &matrix_mode);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity();
  
  gluOrtho2D (-screen.dim[0], screen.dim[0],
	      -screen.dim[1], screen.dim[1]);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix ();
  glLoadIdentity();
  
  glDisable (GL_DEPTH_TEST);
  
  
  for (int l = 0; l < 3; l++)
    {
      int voxel = (vis->screen_voxels * vis->mouse.x[l]) / vis->viewport_pixels[l];
      
      x1[l] = screen.dim[l] * (-1.0 + (voxel<<1) / (double)vis->screen_voxels);
      x2[l] = x1[l] + (screen.dim[l] * (2.0 / vis->screen_voxels));
    }
  glColor3d (0.0, 0.0, 1.0);
  glBegin (GL_LINE_LOOP);
  
  glVertex2d (x1[0], x1[1]);
  glVertex2d (x1[0], x2[1]);
  glVertex2d (x2[0], x2[1]);
  glVertex2d (x2[0], x1[1]);
  
  glEnd ();
  glEnable (GL_DEPTH_TEST);
  
  
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix ();
  
  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
  glMatrixMode (matrix_mode);
}


#ifndef MESH
void visVisualiseSelectedSlice (Vis *vis)
{
  double x, y;
  double dx, dy;
  double grey, delta_grey;
  
  int matrix_mode;
  int voxel[3];
  
  
  glGetIntegerv (GL_MATRIX_MODE, &matrix_mode);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity();
  
  gluOrtho2D (-screen.dim[0], screen.dim[0],
	      -screen.dim[1], screen.dim[1]);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix ();
  glLoadIdentity();
  
  glDisable (GL_DEPTH_TEST);
  glClear (GL_COLOR_BUFFER_BIT);
  
  dx = 2.0*screen.dim[0]/vis->input_voxels[0];
  dy = 2.0*screen.dim[1]/vis->input_voxels[1];
  
  delta_grey = vis->grey_max - vis->grey_min;
  
  voxel[2] = vis->selected_voxel[2];
  
  y = 0.5 * dy - screen.dim[1];
  
  for (int j = 0; j < vis->input_voxels[1] - 1; j++)
    {
      glBegin (GL_QUAD_STRIP);
      
      x = 0.5 * dx - screen.dim[0];
      
      for (int i = 0; i < vis->input_voxels[0]; i++)
	{
	  voxel[0] = i;
	  voxel[1] = j;
	  
	  grey = (vis->voxel[ VoxelId(voxel,vis->input_voxels) ]-vis->grey_min) / delta_grey;
	  
	  glColor3d (grey, grey, grey);	  
  	  glVertex2d (x, y);
  	  
	  voxel[0] = i;
	  voxel[1] = j + 1;
	  
	  grey = (vis->voxel[ VoxelId(voxel,vis->input_voxels) ]-vis->grey_min) / delta_grey;
	  
	  glColor3d (grey, grey, grey);
  	  glVertex2d (x, y + dy);
	  
	  x += dx;
  	}
      y += dy;
      
      glEnd ();
    }
  visVisualiseVisData (vis);
  
  glutSwapBuffers ();
  
  glEnable (GL_DEPTH_TEST);
  
  
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix ();
  
  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
  glMatrixMode (matrix_mode);
}


void visUpdateSegmentation (Vis *vis)
{
  double temp = vis->selected_grey;
  
  if (vis->mouse.dy > 0)
    {
      vis->selected_grey += 0.1 * (vis->mouse.dy * vis->mouse.dy);
    }
  else
    {
      vis->selected_grey -= 0.1 * (vis->mouse.dy * vis->mouse.dy);
    }
  if (vis->coords[C] > 0)
    {
      if (vis->selected_grey < temp)
	{
	  if (segUpdateSegmentation (vis) == SUCCESS)
	    {
	      visCalculateSceneCenter (vis);
	      visProjection (vis);
	    }
	  else
	    {
	      vis->selected_grey = temp;
	      segSegmentation (vis);
	    }
	}
      else if (segSegmentation (vis) == SUCCESS)
	{
	  visCalculateSceneCenter (vis);
	  visProjection (vis);
	}
      else
	{
	  vis->selected_grey = temp;
	  segSegmentation (vis);
	}
    }
  else if (segSegmentation (vis) == SUCCESS)
    {
      visCalculateSceneCenter (vis);
      visProjection (vis);
    }
  else
    {
      vis->selected_grey = temp;
      segSegmentation (vis);
    }
}
#else // MESH


void visVisualiseHitData (Hit *first_hit, Hit *second_hit, Vis *vis)
{
  double x[3];
  
  int dir, l;
  
  Ray ray;
  
  
  for (l = 0; l < 3; l++)
    x[l] = 0.5 * (first_hit->pos[l] + second_hit->pos[l]);
  
  glLineWidth (2.0);
  glBegin (GL_LINES);
  
  glColor3d (1.0, 0.0, 0.0);
  
  for (dir = 0; dir < 14; dir++)
    {
      for (l = 0; l < 3; l++)
	ray.org[l] = x[l];
      
      for (l = 0; l < 3; l++)
	ray.dir[l] = (double)e[ dir*3+l ];
      
      ray.t_max = 1.0e+30;
      ray.t_near = 0.0;
      
      second_hit->previous_triangle_id = -1;
      
      if (rtTraceRay (&ray, second_hit, &vis->mesh) == SUCCESS)
	{
	  glVertex3dv (x);
	  glVertex3dv (second_hit->pos);
	}
      else
	{
	  glVertex3dv (x);
	  glVertex3d (ray.org[0] + ray.t_max * ray.dir[0],
		      ray.org[1] + ray.t_max * ray.dir[1],
		      ray.org[2] + ray.t_max * ray.dir[2]);
	}
    }
  glEnd ();
}


void visVisualiseMesh (Vis *vis)
{
  double dx[3];
  double grey;
  
  int l, m, n;
  
  
  glLineWidth (1.0);
  
  for (n = 0; n < vis->mesh.triangles; n++)
    {
      for (l = 0; l < 3; l++)
	dx[l] = viewpoint.pos[l] - vis->mesh.triangle[n].v[0].pos[l];
      
      grey = fabs(ScalarProd (dx, vis->mesh.triangle[n].nor) / sqrt(ScalarProd (dx, dx)));
      
      glBegin (GL_LINE_LOOP);
      glColor3d (grey, grey, grey);
      
      for (m = 0; m < 3; m++)
	glVertex3dv (vis->mesh.triangle[n].v[m].pos);
      
      glEnd ();
    }
}
#endif // MESH


void visVisualiseFluidSitesWithPoints (Vis *vis)
{
  double seconds = myClock ();
  
  double x1[3], x2[3];
  double col[3];
  double scale;
  
  int voxel[2];
  int l, n;
  
  ScreenVoxel *screen_voxel_p;
  
  
#ifndef MESH
  if (vis->menu.option == CHANGE_THRESHOLD)
    {
      visUpdateSegmentation (vis);
    }
#endif
  glShadeModel (GL_FLAT);
  
  glPointSize (2.0);
  glBegin (GL_POINTS);
  
  for (n = 0; n < vis->screen_voxels * vis->screen_voxels; n++)
    {
      vis->screen_voxel[n].site[0] = -1;
      vis->screen_voxel[n].z[SITE] = 0.0;
    }
#ifndef MESH
  scale = 1.0 / vis->res_factor;
#else
  scale = vis->mesh.voxel_size / vis->res_factor;
#endif
  for (n = 0; n < vis->coords[C]; n++)
    {
#ifndef MESH
      for (l = 0; l < 3; l++)
	x1[l] = (vis->coord[C][n].x[l] + 0.5) * scale - vis->half_dim[l];
#else
      for (l = 0; l < 3; l++)
	x1[l] = vis->seed_pos[l] + (vis->coord[C][n].x[l] - vis->seed_site[l]) * scale;
#endif
      visColorPalette (vis->coord[C][n].iters, vis->res_factor, col);
      
      glColor3dv (col);
      glVertex3dv (x1);
      
      visProject (x1, x2);
      
      voxel[0] = (int)((vis->screen_voxels / (2.0 * screen.dim[0])) * (x2[0] + screen.dim[0]));
      voxel[1] = (int)((vis->screen_voxels / (2.0 * screen.dim[1])) * (x2[1] + screen.dim[1]));
      
      if (voxel[0] < 0 || voxel[0] >= vis->screen_voxels ||
	  voxel[1] < 0 || voxel[1] >= vis->screen_voxels)
	{
	  continue;
	}
      screen_voxel_p = &vis->screen_voxel[ voxel[0]*vis->screen_voxels+voxel[1] ];
      
      if (x2[2] > screen_voxel_p->z[SITE])
	{
	  screen_voxel_p->z[SITE] = x2[2];
	  
	  for (l = 0; l < 3; l++)
	    screen_voxel_p->site[l] = vis->coord[C][n].x[l];
	}
    }
  glEnd ();
  
  glShadeModel (GL_SMOOTH);
  
  vis->fps = 1.0 / (myClock () - seconds);
  
  vis->mode = 2;
}


void visVisualiseFluidSitesWithCubes (Vis *vis)
{
  double seconds = myClock ();
  
  double x[3], col[3];
  double scale;
  
  int l, n;
  
  
  glShadeModel (GL_FLAT);
#ifndef MESH
  scale = 1.0 / vis->res_factor;
#else
  scale = vis->mesh.voxel_size / vis->res_factor;
#endif
  for (n = 0; n < vis->coords[C]; n++)
    {
#ifndef MESH
      for (l = 0; l < 3; l++)
	x[l] = (vis->coord[C][n].x[l] + 0.5) * scale - vis->half_dim[l];
#else
      for (l = 0; l < 3; l++)
	x[l] = vis->seed_pos[l] + (vis->coord[C][n].x[l] - vis->seed_site[l]) * scale;
#endif
      visColorPalette (vis->coord[C][n].iters, vis->res_factor, col);
      
      glColor3dv (col);
      glPushMatrix ();
      glTranslatef (x[0], x[1], x[2]);
      glCallList (1);
      glPopMatrix ();
    }
  glShadeModel (GL_SMOOTH);
  
  vis->fps = 1.0 / (myClock () - seconds);
}


void visVisualiseSystem (Vis *vis)
{
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
#ifdef MESH
      visVisualiseMesh (vis);
#endif
  if (vis->mode == 1)
    {
      visVisualiseFluidSitesWithPoints (vis);
    }
  else
    {
      visVisualiseFluidSitesWithCubes (vis);
    }
  if (vis->mouse.b_id == INLET_BOUNDARY || vis->mouse.b_id == OUTLET_BOUNDARY)
    {
      visVisualiseTrianglePars (vis->mouse.b_id, vis->mouse.t_id, vis);
    }
  else if (vis->mouse.b_id < 0)
    {
#ifndef MESH
      ScreenVoxel *voxel_p = visScreenVoxelPtr (vis->mouse.x, vis);
      
      if (voxel_p->site[0] != -1)
	{
	  visVisualiseSiteData (voxel_p->site, vis);
	}
#else
      visVisualiseSiteData (vis);
#endif
    }
  visVisualiseVisData (vis);
  visVisualiseTriangles (vis);
  visVisualiseDiscs (vis);
  
  glutSwapBuffers ();
}


void GLUTCALLBACK Visualise (void)
{
  void (*VisualisePtr[2]) (Vis *vis);
#ifndef MESH
  VisualisePtr[0] = visVisualiseSelectedSlice;
#endif
  VisualisePtr[1] = visVisualiseSystem;
  
  (*VisualisePtr[ min(1, vis.mode) ]) (&vis);
}


void visInitBoundaries (Vis *vis)
{
  vis->boundary[ INLET_BOUNDARY ].triangle = new BoundaryTriangle [1U<<BOUNDARY_ID_BITS];
  vis->boundary[ INLET_BOUNDARY ].triangles = 0;
  
  vis->boundary[ OUTLET_BOUNDARY ].triangle = new BoundaryTriangle [1U<<BOUNDARY_ID_BITS];
  vis->boundary[ OUTLET_BOUNDARY ].triangles = 0;
  
  vis->boundary[ WALL_BOUNDARY ].triangle = new BoundaryTriangle [1U<<BOUNDARY_ID_BITS];
  vis->boundary[ WALL_BOUNDARY ].triangles = 0;
  
  vis->screen_voxels = 100;
  vis->screen_voxel = new ScreenVoxel[vis->screen_voxels * vis->screen_voxels];
}


void visEndBoundaries (Vis *vis)
{
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      free(vis->boundary[n].triangle);
      vis->boundary[n].triangles = 0;
    }
}


void visInit (int argc, char *argv[], Vis *vis)
{
  int is_checkpoint;
  
  
#ifndef MESH
  if (argc != 8)
    {
      visUsage(argv[0]);
      exit(1);
    }
  vis->input_path    = argv[1];
  vis->output_config = argv[2];
  vis->output_pars   = argv[3];
  vis->checkpoint    = argv[4];
  vis->slice_size    = atof(argv[5]);
  vis->pixel_size    = atof(argv[6]);
  is_checkpoint      = atoi(argv[7]);
  
  vis->mode = 0;
#else
  if (argc != 7)
    {
      visUsage(argv[0]);
      exit(1);
    }
  vis->input_file    = argv[1];
  vis->output_config = argv[2];
  vis->output_pars   = argv[3];
  vis->output_coords = argv[4];
  vis->checkpoint    = argv[5];
  is_checkpoint      = atoi(argv[6]);
  
  vis->mode = 1;
#endif
  
  visInitBoundaries (vis);
  
  vis->res_factor = 1;
  ioReadConfig (vis);
  
  if (is_checkpoint)
    {
      ioReadCheckpoint (vis);
    }
  else
    {
      vis->res_factor = 1;
    }
  
#ifdef MESH
  rtInitRayTracing (&vis->mesh);
#endif
  
  if (!is_checkpoint)
    {
      vis->tot_sites = 0;
      vis->coords[C] = 0;
      vis->segmentation_time = 0.0;
#ifndef MESH
      vis->selected_voxel[2] = 0;
      vis->scene_center[0] = 0.0;
      vis->scene_center[1] = 0.0;
      vis->scene_center[2] = 0.0;
#else
      vis->seed_site[0] = -1;
#endif
    }
  else
    {
      segSegmentation (vis);
#ifndef MESH
      visCalculateSceneCenter (vis);
#endif
    }
  vis->viewport_pixels[0] = 800;
  vis->viewport_pixels[1] = 800;
  
  vis->background[0] = 1.0;
  vis->background[1] = 1.0;
  vis->background[2] = 1.0;
  vis->longitude     = 0.0;
  vis->latitude      = 0.0;
  vis->zoom          = 1.0;
  
  vis->ortho[0] = 0.5 * vis->system_size;
  vis->ortho[1] = 0.5 * vis->system_size;
  vis->viewpoint_radius = 2.0 * vis->system_size;
  vis->viewport_radius = 0.5 * vis->viewpoint_radius;
  
  glutInit (&argc, argv);
  
  visOpenWindow (vis->viewport_pixels[0], vis->viewport_pixels[1]);
  visProjection (vis);
  
  visCreateCubeDisplayList (vis);
  
  menuCreateMenu (vis);
}


void visEnd (Vis *vis)
{
  visDeleteCubeDisplayList ();
  
  free(vis->screen_voxel);
  
  for (int i = 0; i < COORD_BUFFERS; i++)
    free(vis->coord[i]);
  
  free(vis->stack_site);
  
  free(vis->block);
#ifndef MESH
  free(vis->voxel);
#endif
}
