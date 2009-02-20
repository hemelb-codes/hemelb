#include "vis.h"


void visColorPalette (int iters, int res_factor, float col[3])
{
  float t;
  
  
  iters = (iters / res_factor) % 200;
  
  if (iters > 100) iters = 200 - iters;
  
  t = 0.01F * (float)iters;
  
  if (t > 1.F)
    {
      col[0] = 1.F;
      col[1] = 0.F;
      col[2] = 0.F;
    }
  else if (t > 0.0F && t <= 0.25F)
    {
      col[0] = 0.F;
      col[1] = 4.F * t;
      col[2] = 1.F;
    }
  else if (t > 0.25F && t <= 0.5F)
    {
      col[0] = 0.F;
      col[1] = 1.F;
      col[2] = 2.F - 4.F * t;
    }
  else if (t > 0.5F && t <= 0.75F)
    {
      col[0] = 4.F * (t - 0.5F);
      col[1] = 1.F;
      col[2] = 0.F;
    }
  else if (t > 0.75F && t <= 1.0F)
    {
      col[0] = 1.F;
      col[1] = 4.F - 4.F * t;
      col[2] = 0.F;
    }
  else
    {
      col[0] = 0.F;
      col[1] = 0.F;
      col[2] = 1.F;
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
  glutInitWindowSize (pixels_x - 1, pixels_y - 1);
  
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
  float temp;
  
  
  for (int l = 0; l < 3; l++)
    screen.col[l] = vis->background[l];
  
  screen.dim[0] = vis->ortho[0] / vis->zoom;
  screen.dim[1] = vis->ortho[1] / vis->zoom;
  
  screen.pixels[0] = vis->viewport_pixels[0];
  screen.pixels[1] = vis->viewport_pixels[1];
  
  temp = vis->longitude * DEG_TO_RAD;
  
  viewpoint.sin_longitude = sinf(temp);
  viewpoint.cos_longitude = cosf(temp);
  
  temp = vis->latitude * DEG_TO_RAD;
  
  viewpoint.sin_latitude = sinf(temp);
  viewpoint.cos_latitude = cosf(temp);
  
  temp = vis->viewpoint_radius * viewpoint.cos_latitude;
  
  viewpoint.pos[0] = temp * viewpoint.sin_longitude;
  viewpoint.pos[1] = vis->viewpoint_radius * viewpoint.sin_latitude;
  viewpoint.pos[2] = temp * viewpoint.cos_longitude;
  
  for (int l = 0; l < 3; l++)
    viewpoint.pos[l] += vis->scene_center[l];
  
  viewpoint.dist = vis->viewport_radius;
  
  temp = vis->viewport_radius / vis->viewpoint_radius;
  
  for (int l = 0; l < 3; l++)
    screen.ctr[l] = viewpoint.pos[l] + temp*(vis->scene_center[l] - viewpoint.pos[l]);
  
  screen.zoom = vis->zoom;
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  
  glFrustum (-screen.dim[0], screen.dim[0],
	     -screen.dim[1], screen.dim[1],
	     vis->viewport_radius, 1.0e+30F);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  if ((int)((fabsf(vis->latitude) + 90.0) / 180.0) % 2 == 0)
    {
      gluLookAt (viewpoint.pos[0], viewpoint.pos[1], viewpoint.pos[2], 
		 screen.ctr[0], screen.ctr[1], screen.ctr[2],
		 0.0F, 1.0F, 0.0F);
    }
  else
    {
      gluLookAt (viewpoint.pos[0], viewpoint.pos[1], viewpoint.pos[2], 
		 screen.ctr[0], screen.ctr[1], screen.ctr[2],
		 0.0F, -1.0F, 0.0F);
    }
  glClearColor (screen.col[0], screen.col[1], screen.col[2], 0.F);
}


void visProject (float px1[3], float px2[3])
{
  float x1[3], x2[3];
  
  
  for (int l = 0; l < 3; l++)
    x1[l] = px1[l] - viewpoint.pos[l];
  
  AntiRotate (x1,
	      viewpoint.sin_longitude, viewpoint.cos_longitude,
	      viewpoint.sin_latitude, viewpoint.cos_latitude,
	      x2);
  
  float temp = viewpoint.dist / (px2[2] = -x2[2]);
  
  px2[0] = temp * x2[0];
  px2[1] = temp * x2[1];
}


void visAntiProject (float px1[3], float px2[3])
{
  float x1[3], x2[3];
  float temp;
  
  
  temp = px1[2] / viewpoint.dist;
  
  x1[0] = temp * px1[0];
  x1[1] = temp * px1[1];
  x1[2] = -px1[2];
  
  Rotate (x1,
	  viewpoint.sin_longitude, viewpoint.cos_longitude,
	  viewpoint.cos_latitude, viewpoint.cos_latitude,
	  x2);
  
  for (int l = 0; l < 3; l++)
    px2[l] = x2[l] + viewpoint.pos[l];
}


void visRescaleViewpoint (float scale, Vis *vis)
{
  vis->ortho[0] *= scale;
  vis->ortho[1] *= scale;
  vis->viewpoint_radius *= scale;
  vis->viewport_radius = 0.5F * vis->viewpoint_radius;
  
  visProjection (vis);
}


void visCalculateSceneCenter (Vis *vis)
{
  int block_min[3], block_max[3];
  int b[3];
  
  
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
  for (int l = 0; l < 3; l++)
    vis->scene_center[l] = BLOCK_SIZE*0.5F*(block_min[l] + block_max[l]) - vis->half_dim[l];
}


void visCreateCubeDisplayList (void)
{
  glNewList (1, GL_COMPILE);
  
  glBegin (GL_QUAD_STRIP);
  glVertex3f (0.F, 0.F, 0.F);
  glVertex3f (0.F, 1.F, 0.F);
  glVertex3f (0.F, 0.F, 1.F);
  glVertex3f (0.F, 1.F, 1.F);
  glVertex3f (1.F, 0.F, 1.F);
  glVertex3f (1.F, 1.F, 1.F);
  glVertex3f (1.F, 0.F, 0.F);
  glVertex3f (1.F, 1.F, 0.F);
  glVertex3f (0.F, 0.F, 0.F);
  glVertex3f (0.F, 1.F, 0.F);
  glEnd ();
  
  glBegin (GL_QUADS);
  glVertex3f (0.F, 0.F, 0.F);
  glVertex3f (0.F, 0.F, 1.F);
  glVertex3f (1.F, 0.F, 1.F);
  glVertex3f (1.F, 0.F, 0.F);
  glVertex3f (0.F, 1.F, 0.F);
  glVertex3f (0.F, 1.F, 1.F);
  glVertex3f (1.F, 1.F, 1.F);
  glVertex3f (1.F, 1.F, 0.F);
  glEnd ();
  
  glEndList ();
}


void visDeleteCubeDisplayList (void)
{
  glDeleteLists (1, 1);
}


void visVisualiseString (float r, float g, float b, int x, int y, char *string, void *font)
{
  glColor3f (r, g, b);
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
  
  Triangle *t_p;
  
  
  t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  sprintf (pars_string, "pressure = %.1f + %.1f cos(w t + phase)  mmHg, phase = %.1f deg",
	   t_p->pressure_avg, t_p->pressure_amp, t_p->pressure_phs);
  
  visVisualiseString (0.F, 0.F, 0.F, 5, 5,
		      pars_string, GLUT_BITMAP_HELVETICA_12);
}


void visVisualiseSiteData (short int site[3], Vis *vis)
{
  float nor[3];
  float diameter;
  
  char pars_string[256];
  
  
  if (segEstimateNormal (site, nor, vis) == !SUCCESS)
    {
      return;
    }
  segEstimateDiameter (site, nor, &diameter, vis);
  
  if (diameter >= 0.F)
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
  visVisualiseString (0.F, 0.F, 0.F, 5, 5,
		      pars_string, GLUT_BITMAP_HELVETICA_12);
}


void visVisualiseVisData (Vis *vis)
{
  char data_string[256];
  
  
  if (vis->mode == 0)
    {
      sprintf (data_string, "Pixel size (mm) = %.3f  Slice thickness (mm) = %.3f  Threshold = %.1f  ",
	       vis->pixel_size, vis->slice_size, vis->selected_grey);
      
      visVisualiseString (0.F, 1.F, 0.F,
			  5, vis->viewport_pixels[1]-20,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "Resolution enhancement = %i  #sites = %i  #sup. sites = %i  Seg. time = %.3f s",
	       vis->res_factor, vis->tot_sites, vis->coords[C], vis->segmentation_time);
      
      visVisualiseString (0.F, 1.F, 0.F, 5, vis->viewport_pixels[1]-40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
    }
  else
    {
      sprintf (data_string, "longitude = %.1f", vis->longitude);
      
      visVisualiseString (0.F, 0.F, 0.F, vis->viewport_pixels[0]-110, 55,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "latitude = %.1f", vis->latitude);
      
      visVisualiseString (0.F, 0.F, 0.F, vis->viewport_pixels[0]-110, 40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "zoom = %.1f", vis->zoom);
      
      visVisualiseString (0.F, 0.F, 0.F, vis->viewport_pixels[0]-110, 25,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "Pixel size (mm) = %.3f  Slice thickness (mm) = %.3f  Threshold = %.1f  ",
	       vis->pixel_size, vis->slice_size, vis->selected_grey);
      
      visVisualiseString (0.F, 0.F, 0.F,
			  5, vis->viewport_pixels[1]-20,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "Resolution enhancement = %i  #sites = %i  #sup. sites = %i  Seg. time = %.3f s",
	       vis->res_factor, vis->tot_sites, vis->coords[C], vis->segmentation_time);
      
      visVisualiseString (0.F, 0.F, 0.F, 5, vis->viewport_pixels[1]-40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "FPS = %.2f", vis->fps);
      
      visVisualiseString (0.F, 0.F, 0.F,
			  vis->viewport_pixels[0]-110, 5,
			  data_string, GLUT_BITMAP_HELVETICA_12);
    }
}


void visVisualiseTriangles (Vis *vis)
{
  float x1[3], x2[3], nor[3];
  float length;
  
  Triangle *t_p;
  
  
  glLineWidth (3.F);
  
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      if (n == INLET_BOUNDARY)
	{
	  glColor3f (0.5F, 0.5F, 0.5F);
	}
      else if (n == OUTLET_BOUNDARY)
	{
	  glColor3f (0.0F, 0.5F, 0.0F);
	}
      else if (n == WALL_BOUNDARY)
	{
	  glColor3f (0.0F, 0.0F, 0.5F);
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
	    {
	      glVertex3fv (t_p->v[i].pos);
	    }
	  glEnd ();
	  
	  if (n != INLET_BOUNDARY) continue;
	  
	  glBegin (GL_LINES);
	  
	  for (int l = 0; l < 3; l++)
	    {
	      x1[l] = t_p->pos[l];
	      nor[l] = t_p->nor[l];
	    }
	  length = sqrtf(t_p->d.r2);
	  
	  if (t_p->normal_sign == -1)
	    {
	      for (int l = 0; l < 3; l++)
		nor[l] = -nor[l];
	    }
	  for (int l = 0; l < 3; l++)
	    x2[l] = x1[l] + nor[l] * length;
	  
	  glVertex3fv (x1);
	  glVertex3fv (x2);
	  
	  glEnd ();
	}
    }
}


void visVisualiseDiscs (Vis *vis)
{
  float v[2*37];
  float x1[3], x2[3];
  
  Triangle *t_p;
  
  
  for (int i = 0; i <= 35; i++)
    {
      v[ 2*i   ] = sinf(i * 10 * DEG_TO_RAD);
      v[ 2*i+1 ] = cosf(i * 10 * DEG_TO_RAD);
    }
  v[ 2*36   ] = v[ 0 ];
  v[ 2*36+1 ] = v[ 1 ];
  
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      if (n == INLET_BOUNDARY)
	{
	  glColor3f (0.5F, 0.5F, 0.5F);
	}
      else if (n == OUTLET_BOUNDARY)
	{
	  glColor3f (0.0F, 0.5F, 0.0F);
	}
      else if (n == WALL_BOUNDARY)
	{
	  glColor3f (0.0F, 0.0F, 0.5F);
	}
      for (unsigned int m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  glBegin (GL_TRIANGLE_FAN);
	  
	  t_p = &vis->boundary[ n ].triangle[ m ];
	  
	  glVertex3fv (t_p->pos);
	  
	  for (int i = 0; i <= 36; i++)
	    {
	      x1[0] = sqrtf(t_p->d.r2) * v[ 2*i   ];
	      x1[1] = sqrtf(t_p->d.r2) * v[ 2*i+1 ];
	      x1[2] = 0.0F;
	      
	      Rotate (x1,
		      t_p->d.sin_longitude, t_p->d.cos_longitude,
		      t_p->d.sin_latitude,  t_p->d.cos_latitude,
		      x2);
	      
	      for (int l = 0; l < 3; l++)
		x1[l] = x2[l] + t_p->pos[l];
	      
	      glVertex3fv (x1);
	    }
	  glEnd ();
	}
    }
}


void visVisualiseActiveBoundaryVoxel (Vis *vis)
{
  float x1[2], x2[2];
  
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
      
      x1[l] = screen.dim[l] * (-1.F + (voxel<<1) / (float)vis->screen_voxels);
      x2[l] = x1[l] + (screen.dim[l] * (2.0F / vis->screen_voxels));
    }
  glColor3f (0.F, 0.F, 1.F);
  glBegin (GL_LINE_LOOP);
  
  glVertex2f (x1[0], x1[1]);
  glVertex2f (x1[0], x2[1]);
  glVertex2f (x2[0], x2[1]);
  glVertex2f (x2[0], x1[1]);
  
  glEnd ();
  glEnable (GL_DEPTH_TEST);
  
  
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix ();
  
  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
  glMatrixMode (matrix_mode);
}


void visVisualiseSelectedSlice (Vis *vis)
{
  float x, y;
  float dx, dy;
  float grey, delta_grey;
  
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
  
  dx = 2.0F*screen.dim[0]*vis->viewport_pixels[0] / (vis->input_voxels[0]*(vis->viewport_pixels[0]-2));
  dy = 2.0F*screen.dim[1]*vis->viewport_pixels[1] / (vis->input_voxels[1]*(vis->viewport_pixels[1]-2));
  
  delta_grey = vis->grey_max - vis->grey_min;
  
  voxel[2] = vis->selected_voxel[2];
  
  y = -screen.dim[1];
  
  for (int j = 0; j < vis->input_voxels[1] - 1; j++)
    {
      glBegin (GL_QUAD_STRIP);
      
      x = -screen.dim[0];
      
      for (int i = 0; i < vis->input_voxels[0]; i++)
	{
	  voxel[0] = i;
	  voxel[1] = j;
	  
	  grey = (vis->voxel[ VoxelId(voxel,vis->input_voxels) ]-vis->grey_min) / delta_grey;
	  
	  glColor3f (grey, grey, grey);	  
  	  glVertex2f (x, y);
  	  
	  voxel[0] = i;
	  voxel[1] = j + 1;
	  
	  grey = (vis->voxel[ VoxelId(voxel,vis->input_voxels) ]-vis->grey_min) / delta_grey;
	  
	  glColor3f (grey, grey, grey);
  	  glVertex2f (x, y + dy);
	  
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
  float temp = vis->selected_grey;
  
  if (vis->mouse.dy > 0)
    {
      vis->selected_grey += vis->mouse.dy * vis->mouse.dy;
    }
  else
    {
      vis->selected_grey -= vis->mouse.dy * vis->mouse.dy;
    }
  vis->selected_grey = fmaxf(vis->grey_min, fminf(vis->grey_max, vis->selected_grey));
  
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


void visVisualiseFluidSitesWithPoints (Vis *vis)
{
  double seconds = myClock ();
  
  float x1[3], x2[3];
  float col[3];
  
  int voxel[2];
  int l, n;
  
  ScreenVoxel *screen_voxel_p;
  
  
  if (vis->menu.option == CHANGE_THRESHOLD)
    {
      visUpdateSegmentation (vis);
    }
  glShadeModel (GL_FLAT);
  
  glPointSize (2.0F);
  glBegin (GL_POINTS);
  
  for (n = 0; n < vis->screen_voxels * vis->screen_voxels; n++)
    {
      vis->screen_voxel[n].site[0] = -1;
      vis->screen_voxel[n].z[SITE] = 0.0F;
    }
  for (n = 0; n < vis->coords[C]; n++)
    {
      for (l = 0; l < 3; l++)
	x1[l] = (float)vis->coord[C][n].x[l] + (0.5F - vis->half_dim[l]);
      
      visColorPalette (vis->coord[C][n].iters, vis->res_factor, col);
      
      glColor3fv (col);
      glVertex3fv (x1);
      
      visProject (x1, x2);
      
      voxel[0] = (int)((vis->screen_voxels / (2.0F * screen.dim[0])) * (x2[0] + screen.dim[0]));
      voxel[1] = (int)((vis->screen_voxels / (2.0F * screen.dim[1])) * (x2[1] + screen.dim[1]));
      
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
  
  vis->fps = 1.0F / (myClock () - seconds);
  
  vis->mode = 2;
}


void visVisualiseFluidSitesWithCubes (Vis *vis)
{
  double seconds = myClock ();
  
  float x[3];
  float col[3];
  
  int l, n;
  
  
  glShadeModel (GL_FLAT);
  
  for (n = 0; n < vis->coords[C]; n++)
    {
      for (l = 0; l < 3; l++)
	x[l] = (float)vis->coord[C][n].x[l] + (0.5F - vis->half_dim[l]);
      
      visColorPalette (vis->coord[C][n].iters, vis->res_factor, col);
      
      glColor3fv (col);
      glPushMatrix ();
      glTranslatef (x[0], x[1], x[2]);
      glCallList (1);
      glPopMatrix ();
    }
  glShadeModel (GL_SMOOTH);
  
  vis->fps = 1.0F / (myClock () - seconds);
}


void visVisualiseSystem (Vis *vis)
{
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  if (vis->mode == 1)
    {
      visVisualiseFluidSitesWithPoints (vis);
    }
  else
    {
      visVisualiseFluidSitesWithCubes (vis);
    }
  visVisualiseVisData (vis);
  
  if (vis->mouse.b_id == INLET_BOUNDARY || vis->mouse.b_id == OUTLET_BOUNDARY)
    {
      visVisualiseTrianglePars (vis->mouse.b_id, vis->mouse.t_id, vis);
    }
  else if (vis->mouse.b_id < 0)
    {
      ScreenVoxel *voxel_p = visScreenVoxelPtr (vis->mouse.x, vis);
      
      if (voxel_p->site[0] != -1)
	{
	  visVisualiseSiteData (voxel_p->site, vis);
	}
    }
  visVisualiseVisData (vis);
  visVisualiseTriangles (vis);
  visVisualiseDiscs (vis);
  
  glutSwapBuffers ();
}


void GLUTCALLBACK Visualise (void)
{
  void (*VisualisePtr[2]) (Vis *vis);
  
  VisualisePtr[0] = visVisualiseSelectedSlice;
  VisualisePtr[1] = visVisualiseSystem;
  
  (*VisualisePtr[ min(1, vis.mode) ]) (&vis);
}


void visInitBoundaries (Vis *vis)
{
  vis->boundary[ INLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U<<BOUNDARY_ID_BITS));
  vis->boundary[ INLET_BOUNDARY ].triangles = 0;
  
  vis->boundary[ OUTLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U<<BOUNDARY_ID_BITS));
  vis->boundary[ OUTLET_BOUNDARY ].triangles = 0;
  
  vis->boundary[ WALL_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U<<BOUNDARY_ID_BITS));
  vis->boundary[ WALL_BOUNDARY ].triangles = 0;
  
  vis->screen_voxel = (ScreenVoxel *)malloc(sizeof(ScreenVoxel) * vis->screen_voxels * vis->screen_voxels);
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
  
  
  if (argc == 7)
    {
      is_checkpoint = 0;
    }
  else if (argc == 4)
    {
      is_checkpoint = 1;
    }
  else
    {
      visUsage(argv[0]);
      exit(1);
    }
  
  if (!is_checkpoint)
    {
      vis->input_path    = argv[1];
      vis->output_config = argv[2];
      vis->output_pars   = argv[3];
      vis->checkpoint    = argv[4];
      
      vis->slice_size = atof(argv[5]);
      vis->pixel_size = atof(argv[6]);
      
      vis->res_factor = 1;
      
      vis->selected_voxel[2] = 0;
      
      vis->screen_voxels = 100;
      vis->mode = 0;
      
      ioReadConfig (vis);
      visInitBoundaries (vis);
      
      vis->tot_sites = 0;
      vis->coords[C] = 0;
      vis->segmentation_time = 0.F;
      
      vis->scene_center[0] = 0.F;
      vis->scene_center[1] = 0.F;
      vis->scene_center[2] = 0.F;
    }
  else
    {
      vis->checkpoint    = argv[1];
      vis->output_config = argv[2];
      vis->output_pars   = argv[3];
      
      vis->screen_voxels = 100;
      vis->mode = 1;
      
      ioReadCheckpoint (vis);
      segSegmentation (vis);
      visCalculateSceneCenter (vis);
    }
  
  vis->viewport_pixels[0] = 512;
  vis->viewport_pixels[1] = 512;
  
  vis->background[0] = 1.0F;
  vis->background[1] = 1.0F;
  vis->background[2] = 1.0F;
  
  vis->longitude = 45.F;
  vis->latitude  = 45.F;
  
  vis->zoom = 1.0F;
  
  vis->ortho[0] = 0.5F * vis->system_size;
  vis->ortho[1] = 0.5F * vis->system_size;
  vis->viewpoint_radius = 2.F * vis->system_size;
  vis->viewport_radius = 0.5F * vis->viewpoint_radius;
  
  glutInit (&argc, argv);
  
  visOpenWindow (vis->viewport_pixels[0], vis->viewport_pixels[1]);
  visProjection (vis);
  
  visCreateCubeDisplayList ();
  
  menuCreateMenu (vis);
}


void visEnd (Vis *vis)
{
  visDeleteCubeDisplayList ();
  
  free(vis->screen_voxel);
  
  for (int i = 0; i < COORD_BUFFERS; i++)
    {
      free(vis->coord[i]);
    }
  free(vis->stack_site);
  
  free(vis->block);
  
  free(vis->voxel);
}
