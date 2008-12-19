#define GL_GLEXT_PROTOTYPES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include <rpc/xdr.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


#ifndef GLUTCALLBACK
#define GLUTCALLBACK
#endif

#define EPSILON       1.0e-30
#define DEG_TO_RAD    0.01745329

#define SUCCESS                 1
#define COORDS_BUFFERS          4
#define COORDS_BUFFERS_SIZE_A   100000
#define COORDS_BUFFERS_SIZE_B   10000
#define ITERS_MAX               10
#define VOXELS_X                200
#define VOXELS_Y                200

#define NULL_MENU_OPTION               0
#define SEGMENT_1X                (1<<0)
#define SEGMENT_2X                (1<<1)
#define SEGMENT_3X                (1<<2)
#define SEGMENT_4X                (1<<3)
#define SEGMENT_5X                (1<<4)
#define SEGMENT_6X                (1<<5)
#define CHANGE_THRESHOLD          (1<<6)
#define CHANGE_SLICE              (1<<7)
#define ZOOM_SCENE                (1<<8)
#define ROTATE_SCENE              (1<<9)
#define CREATE_INLET             (1<<10)
#define CREATE_OUTLET            (1<<11)
#define CREATE_WALL              (1<<12)
#define ZOOM_BOUNDARY            (1<<13)
#define ROTATE_BOUNDARY          (1<<14)
#define REVERSE_INLET_NORMAL     (1<<15)
#define DELETE_BOUNDARY          (1<<16)
#define CHANGE_VIS_MODE          (1<<17)
#define SAVE_DATA                (1<<18)
#define QUIT                     (1<<19)
#define ACTIVE                         1


#define visVoxelId(i,j,k,vis) (k * vis->input_slices + j) * vis->input_pixels_x + i


unsigned int SOLID_TYPE  = 0U;
unsigned int FLUID_TYPE  = 1U;
unsigned int INLET_TYPE  = 2U;
unsigned int OUTLET_TYPE = 3U;

unsigned int BOUNDARIES              = 3U;
unsigned int INLET_BOUNDARY          = 0U;
unsigned int OUTLET_BOUNDARY         = 1U;
unsigned int WALL_BOUNDARY           = 2U;

unsigned int SITE_TYPE_BITS       = 2U;
unsigned int BOUNDARY_CONFIG_BITS = 14U;
unsigned int BOUNDARY_DIR_BITS    = 4U;
unsigned int BOUNDARY_ID_BITS     = 10U;

unsigned int BOUNDARY_CONFIG_SHIFT = SITE_TYPE_BITS;
unsigned int BOUNDARY_DIR_SHIFT    = BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
unsigned int BOUNDARY_ID_SHIFT     = BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

unsigned int SITE_TYPE_MASK       = ((1U << SITE_TYPE_BITS) - 1U);
unsigned int BOUNDARY_CONFIG_MASK = ((1U << BOUNDARY_CONFIG_BITS) - 1U) << BOUNDARY_CONFIG_SHIFT;
unsigned int BOUNDARY_DIR_MASK    = ((1U << BOUNDARY_DIR_BITS) - 1U)    << BOUNDARY_DIR_SHIFT;
unsigned int BOUNDARY_ID_MASK     = ((1U << BOUNDARY_ID_BITS) - 1U)     << BOUNDARY_ID_SHIFT;

int e[] = {
   1, 0, 0,
  -1, 0, 0,
   0, 1, 0,
   0,-1, 0,
   0, 0, 1,
   0, 0,-1,
   1, 1, 1,
  -1,-1,-1,
   1, 1,-1,
  -1,-1, 1,
   1,-1, 1,
  -1, 1,-1,
   1,-1,-1,
  -1, 1, 1
};

int inv_dir[] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12};


struct Screen
{
  float col_r;
  float col_g;
  float col_b;
  
  float ctr_x;
  float ctr_y;
  float ctr_z;
  
  float max_x;
  float max_y;
  
  float zoom;
  
  int   pix_x;
  int   pix_y;
};


struct Viewpoint
{
  float pos_x;
  float pos_y;
  float pos_z;
  
  float sin_1;
  float cos_1;
  
  float sin_2;
  float cos_2;
  
  float dist;
};


struct Site
{
  unsigned int data;
  unsigned short int iters;
  
  char label;
  char x, y, z;
};


struct Block
{
  Site *site;
  
  short int is_void;
  short int stack_t_id;
  short int x, y, z;
};


struct Vertex
{
  float pos_x, pos_y, pos_z;
};


struct Disc
{
  float sin_longitude, cos_longitude;
  float sin_latitude, cos_latitude;
  
  float r2;
};


struct Triangle
{
  Vertex v[3];
  
  Disc d;
  
  float pos_x, pos_y, pos_z;
  float nor_x, nor_y, nor_z;
  
  float pressure_avg, pressure_amp, pressure_phs;
  
  int normal_sign;
};


struct Boundary
{
  Triangle *triangle;
  
  int triangles;
};


struct StackTriangle
{
  short int b_id, t_id;
  
  int next;
};


struct Coords
{
  short int x, y, z;
};


struct ScreenVoxel
{
  float v_z;
  float site_z;
  
  short int t_id;
  short int site_i, site_j, site_k;
  
  char b_id;
  char v_id;
};


struct Mouse
{
  short int t_id;
  char b_id;
  char v_id;
  
  short int x, y;
  
  short int state;
};


struct Menu
{
  int id;
  int option;
};


struct Vis
{
  float scale_x, scale_y, scale_z;
  float scale_inv_x, scale_inv_y, scale_inv_z;
  
  int input_pixels_x, input_pixels_y;
  int output_image_pix_x, output_image_pix_y;
  int input_slices, output_slices;
  int pixel_depth;
  
  int sites_x, sites_y, sites_z;
  int blocks_x, blocks_y, blocks_z;
  int blocks;
  int block_size;
  int shift;
  int sites_in_a_block;
  
  int sites;
  int boundary_sites;
  int stack_sites, stack_sites_max;
  int stack_triangles, stack_triangles_max;
  int sites_a;
  int iters;

  
  int selected_pixel_x, selected_pixel_y, selected_slice;
  
  float selected_gray;
  float gray_min, gray_max;
  
  int viewport_pixels_x;
  int viewport_pixels_y;
  
  float background_r;
  float background_g;
  float background_b;
  
  float ortho_x, ortho_y;
  float longitude, latitude;
  float viewpoint_radius;
  float viewport_radius;
  float zoom;
  float scene_center_x;
  float scene_center_y;
  float scene_center_z;
  
  float system_size;
  float dim_x, dim_y, dim_z;
  float half_dim_x, half_dim_y, half_dim_z;
  
  float slice_size, pixel_size;
  float res_factor_par1, res_factor_par2;
  
  int res_factor;
  int screen_voxels;
  int smoothing_range;
  int voxel_di, voxel_dj, voxel_dk;
  int point_size;
  int visualise_boundaries;
  int mode;
  
  float segmentation_time, fps;
  
  
  Site *stack_site;
  
  Block *block;
  
  StackTriangle *stack_triangle;
  
  unsigned short int *medical_data;
  
  Coords *site_coords[COORDS_BUFFERS];
  
  Boundary boundary[4];
  
  ScreenVoxel *screen_voxel;
  
  Mouse mouse;
  
  Menu menu;
  
  int *stored_block;
  
  char *input_path;
  char *output_config;
  char *output_pars;
  char *checkpoint;
  
  vector<string> file_list;
};


Screen screen;

Viewpoint viewpoint;

Vis vis;


void visInit (int argc, char *argv[], Vis *vis);
void visEnd (Vis *vis);
int visEstimateBoundaryNormal (int site_i, int site_j, int site_k, float nor[], Vis *vis);
void visEstimateDiameter (int site_i, int site_j, int site_k, float nor[], float *diameter, Vis *vis);


int min (int a, int b)
{
  if (a < b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


int max (int a, int b)
{
  if (a > b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


int nint (float a)
{
  if (a > (int)(a + 0.5F))
    {
      return 1 + (int)a;
    }
  else
    {
      return (int)a;
    }
}


double myClock ()
{
  return (double)clock () * (1. / (double)CLOCKS_PER_SEC);
}


void visRotate (float x0, float y0, float z0,
		float longitude, float latitude,
		float *x1, float *y1, float *z1)
{
  float sn1, cs1;
  float sn2, cs2;
  float temp;
  
  
  sn1 = sinf(longitude);
  cs1 = cosf(longitude);
  sn2 = sinf(latitude);
  cs2 = cosf(latitude);

  temp = cs2 * z0 - sn2 * y0;

  *x1 = sn1 * temp + cs1 * x0;
  *y1 = sn2 * z0   + cs2 * y0;
  *z1 = cs1 * temp - sn1 * x0;
}


void visRotate (float x0, float y0, float z0,
		float sn1, float cs1, float sn2, float cs2,
		float *x1, float *y1, float *z1)
{
  float temp;
  
  
  temp = cs2 * z0 - sn2 * y0;

  *x1 = sn1 * temp + cs1 * x0;
  *y1 = sn2 * z0   + cs2 * y0;
  *z1 = cs1 * temp - sn1 * x0;
}


void visAntiRotate (float x0, float y0, float z0,
		    float longitude, float latitude,
		    float *x1, float *y1, float *z1)
{
  float sn1, cs1;
  float sn2, cs2;
  float temp;
  
  
  sn1 = sinf(longitude);
  cs1 = cosf(longitude);
  sn2 = sinf(latitude);
  cs2 = cosf(latitude);
  
  temp = cs1 * z0 + sn1 * x0;

  *x1 = cs1 * x0   - sn1 * z0;
  *y1 = cs2 * y0   - sn2 * temp;
  *z1 = cs2 * temp + sn2 * y0;
}


void visAntiRotate (float x0, float y0, float z0,
		    float sn1, float cs1, float sn2, float cs2,
		    float *x1, float *y1, float *z1)
{
  float temp;
  
  
  temp = cs1 * z0 + sn1 * x0;
  
  *x1 = cs1 * x0   - sn1 * z0;
  *y1 = cs2 * y0   - sn2 * temp;
  *z1 = cs2 * temp + sn2 * y0;
}


void visTriangleCenter (Triangle *t_p, float *cx, float *cy, float *cz)
{
  *cx = (1.F/3.F) * (t_p->v[0].pos_x + t_p->v[1].pos_x + t_p->v[2].pos_x);
  *cy = (1.F/3.F) * (t_p->v[0].pos_y + t_p->v[1].pos_y + t_p->v[2].pos_y);
  *cz = (1.F/3.F) * (t_p->v[0].pos_z + t_p->v[1].pos_z + t_p->v[2].pos_z);
}


void visTriangleNormal (Triangle *t_p, float *nx, float *ny, float *nz)
{
  float dx1, dy1, dz1;
  float dx2, dy2, dz2;
  
  float temp;
  
  
  dx1 = t_p->v[1].pos_x - t_p->v[0].pos_x;
  dy1 = t_p->v[1].pos_y - t_p->v[0].pos_y;
  dz1 = t_p->v[1].pos_z - t_p->v[0].pos_z;
  dx2 = t_p->v[2].pos_x - t_p->v[0].pos_x;
  dy2 = t_p->v[2].pos_y - t_p->v[0].pos_y;
  dz2 = t_p->v[2].pos_z - t_p->v[0].pos_z;
  
  *nx = dy2 * dz1 - dz2 * dy1;
  *ny = dz2 * dx1 - dx2 * dz1;
  *nz = dx2 * dy1 - dy2 * dx1;
  
  temp = 1.F / fmaxf(1.0e-30F, sqrtf(*nx * *nx + *ny * *ny + *nz * *nz));
  
  *nx *= temp;
  *ny *= temp;
  *nz *= temp;
}


float visTriangleRadius (Triangle *t_p)
{
  float cx, cy, cz;
  float dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3;
  float mx1, my1, mz1, mx2, my2, mz2, mx3, my3, mz3;
  float d1, d2, d3;
  
  
  mx1 = 0.5F * (t_p->v[1].pos_x + t_p->v[0].pos_x);
  my1 = 0.5F * (t_p->v[1].pos_y + t_p->v[0].pos_y);
  mz1 = 0.5F * (t_p->v[1].pos_z + t_p->v[0].pos_z);
  mx2 = 0.5F * (t_p->v[2].pos_x + t_p->v[0].pos_x);
  my2 = 0.5F * (t_p->v[2].pos_y + t_p->v[0].pos_y);
  mz2 = 0.5F * (t_p->v[2].pos_z + t_p->v[0].pos_z);
  mx3 = 0.5F * (t_p->v[2].pos_x + t_p->v[1].pos_x);
  my3 = 0.5F * (t_p->v[2].pos_y + t_p->v[1].pos_y);
  mz3 = 0.5F * (t_p->v[2].pos_z + t_p->v[1].pos_z);
  
  visTriangleCenter (t_p, &cx, &cy, &cz);
  
  dx1 = mx1 - cx;
  dy1 = my1 - cy;
  dz1 = mz1 - cz;
  dx2 = mx2 - cx;
  dy2 = my2 - cy;
  dz2 = mz2 - cz;
  dx3 = mx3 - cx;
  dy3 = my3 - cy;
  dz3 = mz3 - cz;
  
  d1 = dx1 * dx1 + dy1 * dy1 + dy1 * dy1;
  d2 = dx2 * dx2 + dy2 * dy2 + dy2 * dy2;
  d3 = dx3 * dx3 + dy3 * dy3 + dy3 * dy3;
  
  
  if (d1 < d2 && d1 < d3)
    {
      return sqrtf(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
    }
  else if (d2 < d1 && d2 < d3)
    {
      return sqrtf(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
    }
  else
    {
      return sqrtf(dx3 * dx3 + dy3 * dy3 + dz3 * dz3);
    }
}


float visTriangleArea (Triangle *t_p)
{
  float dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3;
  float a, b, c, p;
  
  
  dx1 = t_p->v[1].pos_x - t_p->v[0].pos_x;
  dy1 = t_p->v[1].pos_y - t_p->v[0].pos_y;
  dz1 = t_p->v[1].pos_z - t_p->v[0].pos_z;
  dx2 = t_p->v[2].pos_x - t_p->v[0].pos_x;
  dy2 = t_p->v[2].pos_y - t_p->v[0].pos_y;
  dz2 = t_p->v[2].pos_z - t_p->v[0].pos_z;
  dx3 = t_p->v[2].pos_x - t_p->v[1].pos_x;
  dy3 = t_p->v[2].pos_y - t_p->v[1].pos_y;
  dz3 = t_p->v[2].pos_z - t_p->v[1].pos_z;
  
  a = sqrtf(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  b = sqrtf(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
  c = sqrtf(dx3 * dx3 + dy3 * dy3 + dz3 * dz3);
  
  p = 0.5F * (a + b + c);
  
  return sqrtf(p * (p - a) * (p - b) * (p - c));
}


int visRayVsTriangle (Triangle *t_p,
		      float px, float py, float pz,
		      float nx, float ny, float nz,
		      float t_max, float *t)
{
  float ex1, ey1, ez1, ex2, ey2, ez2;
  float det;
  float x1, y1, z1, x2, y2, z2, x3, y3, z3;
  float v, w;
  
  
  ex1 = t_p->v[1].pos_x - t_p->v[0].pos_x;
  ey1 = t_p->v[1].pos_y - t_p->v[0].pos_y;
  ez1 = t_p->v[1].pos_z - t_p->v[0].pos_z;
  ex2 = t_p->v[2].pos_x - t_p->v[0].pos_x;
  ey2 = t_p->v[2].pos_y - t_p->v[0].pos_y;
  ez2 = t_p->v[2].pos_z - t_p->v[0].pos_z;
  
  x1 = ny * ez2 - nz * ey2;
  y1 = nz * ex2 - nx * ez2;
  z1 = nx * ey2 - ny * ex2;
  
  det = ex1 * x1 + ey1 * y1 + ez1 * z1;
  
  if (det > -1.e-3F && det < 1.e-30F)
    {
      return !SUCCESS;
    }
  det = 1.F / det;
  
  x2 = px - t_p->v[0].pos_x;
  y2 = py - t_p->v[0].pos_y;
  z2 = pz - t_p->v[0].pos_z;
  
  v = (x2 * x1 + y2 * y1 + z2 * z1) * det;
  
  if (v < 0.F || v > 1.F)
    {
      return !SUCCESS;
    }
  
  x3 = y2 * ez1 - z2 * ey1;
  y3 = z2 * ex1 - x2 * ez1;
  z3 = x2 * ey1 - y2 * ex1;
  
  w = (nx * x3 + ny * y3 + nz * z3) * det;
  
  if (w < 0.F || v + w > 1.F)
    {
      return !SUCCESS;
    }

  if ((*t = (ex2 * x3 + ey2 * y3 + ez2 * z3) * det) < EPSILON)
    {
      return !SUCCESS;
    }
  if (*t > t_max)
    {
      return !SUCCESS;
    }
  return SUCCESS;
}


int visRayVsDisc (Triangle *t_p,
		  float x1, float y1, float z1,
		  float x2, float y2, float z2,
		  float t_max, float *t)
{
  float x3, y3, z3;
  float x4, y4, z4;
  float x, y;
  
  
  x1 -= t_p->pos_x;
  y1 -= t_p->pos_y;
  z1 -= t_p->pos_z;
  
  x2 -= t_p->pos_x;
  y2 -= t_p->pos_y;
  z2 -= t_p->pos_z;
  
  visAntiRotate (x1, y1, z1,
		 t_p->d.sin_longitude, t_p->d.cos_longitude,
		 t_p->d.sin_latitude,  t_p->d.cos_latitude,
		 &x3, &y3, &z3);
  
  visAntiRotate (x2, y2, z2,
		 t_p->d.sin_longitude, t_p->d.cos_longitude,
		 t_p->d.sin_latitude,  t_p->d.cos_latitude,
		 &x4, &y4, &z4);
  
  if ((z3 > 0.F && z4 > 0.F) ||
      (z3 < 0.F && z4 < 0.F))
    {
      return !SUCCESS;
    }
  *t = z3 / (z3 - z4);
  
  if (*t >= t_max)
    {
      return !SUCCESS;
    }
  x = (1.F - *t) * x3 + *t * x4;
  y = (1.F - *t) * y3 + *t * y4;
  
  if (x * x + y * y > t_p->d.r2)
    {
      return !SUCCESS;
    }
  return SUCCESS;
}


void visOpenWindow (int pixels_x, int pixels_y)
{
  screen.pix_x = pixels_x;
  screen.pix_y = pixels_y;
  
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
  
  
  screen.col_r = vis->background_r;
  screen.col_g = vis->background_g;
  screen.col_b = vis->background_b;
  
  screen.max_x = vis->ortho_x / vis->zoom;
  screen.max_y = vis->ortho_y / vis->zoom;
  
  screen.pix_x = vis->viewport_pixels_x;
  screen.pix_y = vis->viewport_pixels_y;
  
  temp = vis->longitude * DEG_TO_RAD;
  
  viewpoint.sin_1 = sinf(temp);
  viewpoint.cos_1 = cosf(temp);
  
  temp = vis->latitude * DEG_TO_RAD;
  
  viewpoint.sin_2 = sinf(temp);
  viewpoint.cos_2 = cosf(temp);
  
  temp = vis->viewpoint_radius * viewpoint.cos_2;
  
  viewpoint.pos_x = temp * viewpoint.sin_1 + vis->scene_center_x;
  viewpoint.pos_y = vis->viewpoint_radius * viewpoint.sin_2 + vis->scene_center_y;
  viewpoint.pos_z = temp * viewpoint.cos_1 + vis->scene_center_z;
  
  viewpoint.dist = vis->viewport_radius;
  
  temp = vis->viewport_radius / vis->viewpoint_radius;
  
  screen.ctr_x = viewpoint.pos_x + temp * (vis->scene_center_x - viewpoint.pos_x);
  screen.ctr_y = viewpoint.pos_y + temp * (vis->scene_center_y - viewpoint.pos_y);
  screen.ctr_z = viewpoint.pos_z + temp * (vis->scene_center_z - viewpoint.pos_z);
  
  screen.zoom = vis->zoom;
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  
  glFrustum (-screen.max_x, screen.max_x,
	     -screen.max_y, screen.max_y,
	     vis->viewport_radius, 1.0e+30F);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  gluLookAt (viewpoint.pos_x, viewpoint.pos_y, viewpoint.pos_z, 
	     screen.ctr_x, screen.ctr_y, screen.ctr_z,
	     0.0F, 1.0F, 0.0F);
  
  glClearColor (screen.col_r, screen.col_g, screen.col_b, 0.F);
}


void visProject (float  px1, float  py1, float  pz1,
		 float *px2, float *py2, float *pz2)
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
  
  temp = viewpoint.dist / (*pz2 = -z2);
  
  *px2 = temp * x2;
  *py2 = temp * y2;
}


void visAntiProject (float  px1, float  py1, float  pz1,
		     float *px2, float *py2, float *pz2)
{
  float x1, y1, z1, x2, y2, z2;
  float temp;
  
  
  temp = pz1 / viewpoint.dist;
  
  x1 = temp * px1;
  y1 = temp * py1;
  z1 = -pz1;
  
  temp = viewpoint.cos_2 * z1 - viewpoint.sin_2 * y1;
  
  x2 = viewpoint.sin_1 * temp + viewpoint.cos_1 * x1;
  y2 = viewpoint.sin_2 * z1   + viewpoint.cos_2 * y1;
  z2 = viewpoint.cos_1 * temp - viewpoint.sin_1 * x1;
  
  *px2 = x2 + viewpoint.pos_x;
  *py2 = y2 + viewpoint.pos_y;
  *pz2 = z2 + viewpoint.pos_z;
}


void visSaveWindowImage (char *file_name)
{
  FILE *ppm_image_file = fopen (file_name, "wb");
  
  int pix_x, pix_y;
  int i, j;
  
  unsigned char *image_data = NULL;  
  unsigned char *image_data_p = NULL;
  unsigned char *row_data;
  
  
  glReadBuffer (GL_FRONT);
  
  pix_x = screen.pix_x;
  pix_y = screen.pix_y;
  
  image_data = (unsigned char *)malloc(sizeof(unsigned char) * pix_x * pix_y * 3);
  
  row_data = (unsigned char *)malloc(sizeof(unsigned char) * pix_x * 3);
  
  image_data_p = image_data;

  for (j = 0; j < pix_y; j++)
    {
      glReadPixels (0, j, pix_x, 1, GL_RGB, GL_UNSIGNED_BYTE, row_data);
      
      for (i = 0; i < pix_x; i++)
	{
	  *image_data_p++ = row_data[ i*3   ];
	  *image_data_p++ = row_data[ i*3+1 ];
	  *image_data_p++ = row_data[ i*3+2 ];
	}
    }
  
  free((unsigned char *)row_data);
  
  fprintf (ppm_image_file, "P6\n%i %i\n255\n", pix_x, pix_y);
  
  for (j = pix_y - 1; j >= 0; j--)
    {
      fwrite (image_data + j * pix_x * 3, 1, pix_x * 3, ppm_image_file);
    }
  
  free((unsigned char *)image_data);
  
  fclose (ppm_image_file);
}


void visColorPalette (int iters, int res_factor, float *r, float *g, float *b)
{
  float t;
  
  
  iters = (iters / res_factor) % 200;
  
  if (iters > 100) iters = 200 - iters;
  
  t = 0.01F * (float)iters;
  
  if (t > 1.F)
    {
      *r = 1.F;
      *g = 0.F;
      *b = 0.F;
    }
  else if (t > 0.0F && t <= 0.25F)
    {
      *r = 0.F;
      *g = 4.F * t;
      *b = 1.F;
    }
  else if (t > 0.25F && t <= 0.5F)
    {
      *r = 0.F;
      *g = 1.F;
      *b = 2.F - 4.F * t;
    }
  else if (t > 0.5F && t <= 0.75F)
    {
      *r = 4.F * (t - 0.5F);
      *g = 1.F;
      *b = 0.F;
    }
  else if (t > 0.75F && t <= 1.0F)
    {
      *r = 1.F;
      *g = 4.F - 4.F * t;
      *b = 0.F;
    }
  else
    {
      *r = 0.F;
      *g = 0.F;
      *b = 1.F;
    }
}


void visFromVoxelToSiteCoords (int pixel_i, int pixel_j, int slice_id,
			       int *site_i, int *site_j, int *site_k, Vis *vis)
{
  *site_i = (int)(pixel_i  * vis->scale_x);
  *site_j = (int)(pixel_j  * vis->scale_y);
  *site_k = (int)(slice_id * vis->scale_z);
}


void visFromSiteToVoxelCoords (int site_i, int site_j, int site_k,
			       int *pixel_i, int *pixel_j, int *slice_id, Vis *vis)
{
  *pixel_i  = (int)(site_i * vis->scale_inv_x);
  *pixel_j  = (int)(site_j * vis->scale_inv_y);
  *slice_id = (int)(site_k * vis->scale_inv_z);
}


void visFromSiteToGridCoords (int site_i, int site_j, int site_k,
			      int *si, int *sj, int *sk,
			      int *block_id, int *site_id, Vis *vis)
{
  int i, j, k;
  
  
  i = site_i >> vis->shift;
  j = site_j >> vis->shift;
  k = site_k >> vis->shift;
  *block_id = (i * vis->blocks_y + j) * vis->blocks_z + k;
  
  *si = site_i - (i << vis->shift);
  *sj = site_j - (j << vis->shift);
  *sk = site_k - (k << vis->shift);
  *site_id = (((*si << vis->shift) + *sj) << vis->shift) + *sk;
}


float visInterpolatedGray (int site_i, int site_j, int site_k, Vis *vis)
{
  float gray[2][2][2];
  float x, y, z;
  
  int voxel_i[2], voxel_j[2], voxel_k[2];
  int i, j, k;
  
  
  voxel_i[0] = (int)(x = (float)site_i * vis->scale_inv_x);
  voxel_j[0] = (int)(y = (float)site_j * vis->scale_inv_y);
  voxel_k[0] = (int)(z = (float)site_k * vis->scale_inv_z);
  
  x -= (float)voxel_i[0];
  y -= (float)voxel_j[0];
  z -= (float)voxel_k[0];
  
  voxel_i[1] = min(voxel_i[0] + 1, vis->input_pixels_x - 1);
  voxel_j[1] = min(voxel_j[0] + 1, vis->input_pixels_y - 1);
  voxel_k[1] = min(voxel_k[0] + 1, vis->input_slices - 1);
  
  for (k = 0; k < 2; k++)
    {
      for (j = 0; j < 2; j++)
        {
          for (i = 0; i < 2; i++)
            {
	      gray[k][j][i] = (float)vis->medical_data[ visVoxelId(voxel_i[i],voxel_j[j],voxel_k[k],vis) ];
	    }
	}
    }
  float gray_00x = (1.F - x) * gray[0][0][0] + x * gray[0][0][1];
  float gray_01x = (1.F - x) * gray[0][1][0] + x * gray[0][1][1];
  float gray_10x = (1.F - x) * gray[1][0][0] + x * gray[1][0][1];
  float gray_11x = (1.F - x) * gray[1][1][0] + x * gray[1][1][1];
  
  float gray_0y = (1.F - y) * gray_00x + y * gray_01x;
  float gray_1y = (1.F - y) * gray_10x + y * gray_11x;
  
  return (1.F - z) * gray_0y + z * gray_1y;
}


//float visInterpolatedGray (int site_i, int site_j, int site_k, Vis *vis)
//{
//  float x, y, z;
//  float dx, dy, dz;
//  float weight, sum1, sum2;
//  float interpolated_gray;
//  
//  int voxel_i, voxel_j, voxel_k;
//  int i, j, k;
//  
//  
//  x = site_i;
//  y = site_j;
//  z = site_k;
//  
//  voxel_i = (int)((float)site_i * vis->scale_inv_x);
//  voxel_j = (int)((float)site_j * vis->scale_inv_y);
//  voxel_k = (int)((float)site_k * vis->scale_inv_z);
//  
//  sum1 = sum2 = 0.F;
//  
//  for (i = max(0, voxel_i-vis->voxel_di); i <= min(voxel_i+vis->voxel_di, vis->input_pixels_x-1); i++)
//    {
//      dx = (float)i * vis->scale_x - x;
//      dx *= dx;
//      
//      for (j = max(0, voxel_j-vis->voxel_dj); j <= min(voxel_j+vis->voxel_dj, vis->input_pixels_y-1); j++)
//	{
//	  dy = (float)j * vis->scale_y - y;
//	  dy *= dy;
//	  
//	  for (k = max(0, voxel_k-vis->voxel_dk); k <= min(voxel_k+vis->voxel_dk, vis->input_slices-1); k++)
//	    {
//	      dz = (float)k * vis->scale_z - z;
//	      dz *= dz;
//	      
//	      weight = dx + dy + dz;
//	      
//	      if (weight - 1.e-6 > vis->res_factor_par1) continue;
//	      
//	      weight = 1.F - weight * vis->res_factor_par2;
//	      sum1 += weight * (float)vis->medical_data[ visVoxelId(i, j, k, vis) ];
//	      sum2 += weight;
//	    }
//	}
//    }
//  interpolated_gray = sum1 / sum2;
//  
//  return interpolated_gray;
//}


void visSetSmoothingRange (int s_range, int res_factor, Vis *vis)
{
  vis->smoothing_range = s_range;
  vis->res_factor = res_factor;
  
  vis->voxel_di = vis->smoothing_range;
  vis->voxel_dj = vis->voxel_di;
  vis->voxel_dk = nint(vis->voxel_di * vis->pixel_size / vis->slice_size);
  
  vis->res_factor_par1 = vis->voxel_di * vis->voxel_di * vis->res_factor * vis->res_factor;
  vis->res_factor_par2 = 1.F / vis->res_factor_par1;
}


Site *visSitePointer (int site_i, int site_j, int site_k, Vis *vis)
{
  int i = site_i >> vis->shift;
  int j = site_j >> vis->shift;
  int k = site_k >> vis->shift;
  
  Block *block_p = &vis->block[(i * vis->blocks_y + j) * vis->blocks_z + k];
  
  if (block_p->site == NULL)
    {
      return NULL;
    }
  else
    {
      int ii = site_i - (i << vis->shift);
      int jj = site_j - (j << vis->shift);
      int kk = site_k - (k << vis->shift);
      
      return &block_p->site[(((ii << vis->shift) + jj) << vis->shift) + kk];
    }
}


Site *visStackSitePointer (Vis *vis)
{
  vis->stack_sites += vis->sites_in_a_block;
  
  if (vis->stack_sites > vis->stack_sites_max)
    {
      vis->stack_sites -= vis->sites_in_a_block;
      return NULL;
    }
  return &vis->stack_site[ vis->stack_sites - vis->sites_in_a_block ];
}


void visUpdateTrianglesGrid (Vis *vis)
{
  float block_size_inv;
  
  int m, n;
  int ii1, jj1, kk1;
  int ii2, jj2, kk2;
  int ii3, jj3, kk3;
  int iia, jja, kka;
  int iib, jjb, kkb;
  int ii0, jj0, kk0;
  
  Triangle *t_p;
  
  Block *block_p;
  
  
  vis->stack_triangles = 0;
  
  for (n = 0; n < vis->blocks; n++)
    {
      vis->block[ n ].stack_t_id = -1;
    }
  block_size_inv = vis->blocks_x / vis->dim_x;
  
  for (n = 0; n < BOUNDARIES; n++)
    {
      for (m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  t_p = &vis->boundary[ n ].triangle[ m ];
	  
	  ii1 = (int)(block_size_inv * (t_p->v[0].pos_x + vis->half_dim_x));
	  jj1 = (int)(block_size_inv * (t_p->v[0].pos_y + vis->half_dim_y));
	  kk1 = (int)(block_size_inv * (t_p->v[0].pos_z + vis->half_dim_z));
	  
	  ii2 = (int)(block_size_inv * (t_p->v[1].pos_x + vis->half_dim_x));
	  jj2 = (int)(block_size_inv * (t_p->v[1].pos_y + vis->half_dim_y));
	  kk2 = (int)(block_size_inv * (t_p->v[1].pos_z + vis->half_dim_z));
	  
	  ii3 = (int)(block_size_inv * (t_p->v[2].pos_x + vis->half_dim_x));
	  jj3 = (int)(block_size_inv * (t_p->v[2].pos_y + vis->half_dim_y));
	  kk3 = (int)(block_size_inv * (t_p->v[2].pos_z + vis->half_dim_z));
	  
	  iia = min(ii1, min(ii2, ii3));
	  iib = max(ii1, max(ii2, ii3));
	  jja = min(jj1, min(jj2, jj3));
	  jjb = max(jj1, max(jj2, jj3));
	  kka = min(kk1, min(kk2, kk3));
	  kkb = max(kk1, max(kk2, kk3));
	  
	  for (ii0 = iia; ii0 <= iib; ii0++)
	    {
	      if (ii0 < 0 || ii0 >= vis->blocks_x) continue;
	      
	      for (jj0 = jja; jj0 <= jjb; jj0++)
		{
		  if (jj0 < 0 || jj0 >= vis->blocks_y) continue;
		  
		  for (kk0 = kka; kk0 <= kkb; kk0++)
		    {
		      if (kk0 < 0 || kk0 >= vis->blocks_z) continue;
		      
		      block_p = &vis->block[ (ii0 * vis->blocks_y + jj0) * vis->blocks_z + kk0 ];
		      
		      if (vis->stack_triangles == vis->stack_triangles_max)
			{
			  vis->stack_triangles_max *= 2;
			  vis->stack_triangle = (StackTriangle *)realloc(vis->stack_triangle,
									 sizeof(StackTriangle) * vis->stack_triangles_max);
			}
		      vis->stack_triangle[ vis->stack_triangles ].b_id = n;
		      vis->stack_triangle[ vis->stack_triangles ].t_id = m;
		      vis->stack_triangle[ vis->stack_triangles ].next = block_p->stack_t_id;
		      block_p->stack_t_id = vis->stack_triangles;
		      ++vis->stack_triangles;
		    }
		}
	    }
	}
    }
}


int visIsSegmentIntercepted (int site_i, int site_j, int site_k, int l,
			     Site *site_p, Vis *vis)
{
  float x1, y1, z1;
  float x2, y2, z2;
  float nx, ny, nz;
  float t, t_max;
  
  int bi, bj, bk;
  int b_id, t_id;
  int stack_t_id;
  int m;
  
  unsigned int n;
  
  
  site_i = (int)((x1 = (float)site_i) + 0.5F * (nx = (float)e[ l*3+0 ]));
  site_j = (int)((y1 = (float)site_j) + 0.5F * (ny = (float)e[ l*3+1 ]));
  site_k = (int)((z1 = (float)site_k) + 0.5F * (nz = (float)e[ l*3+2 ]));
  
  bi = site_i >> vis->shift;
  bj = site_j >> vis->shift;
  bk = site_k >> vis->shift;
  
  stack_t_id = vis->block[ (bi*vis->blocks_y+bj)*vis->blocks_z+bk ].stack_t_id;
  
  if (stack_t_id == -1) return !SUCCESS;
  
  t_max = 1.F;
  b_id = -1;
  
  x2 = (x1 -= vis->half_dim_x) + nx;
  y2 = (y1 -= vis->half_dim_y) + ny;
  z2 = (z1 -= vis->half_dim_z) + nz;
  
  while (stack_t_id != -1)
    {
      n = vis->stack_triangle[ stack_t_id ].b_id;
      m = vis->stack_triangle[ stack_t_id ].t_id;
      
      stack_t_id = vis->stack_triangle[ stack_t_id ].next;
      
      if (!visRayVsDisc (&vis->boundary[ n ].triangle[ m ],
			 x1, y1, z1, x2, y2, z2,
			 t_max, &t))
	{
	  continue;
	}
      t_max = t;
      b_id = n;
      t_id = m;
    }
  if (b_id == -1)
    {
      return !SUCCESS;
    }
  else
    {
      if (b_id == INLET_BOUNDARY)
	{
	  site_p->data = INLET_TYPE | (t_id << BOUNDARY_ID_SHIFT);
	}
      else if (b_id == OUTLET_BOUNDARY)
	{
	  site_p->data = OUTLET_TYPE | (t_id << BOUNDARY_ID_SHIFT);
	}
      return SUCCESS;
    }
}


void visCalculateTriangleData (Triangle *t_p)
{
  float longitude, latitude;
  float r;
  
  
  visTriangleCenter (t_p, &t_p->pos_x, &t_p->pos_y, &t_p->pos_z);
  visTriangleNormal (t_p, &t_p->nor_x, &t_p->nor_y, &t_p->nor_z);
  
  longitude = atan2f(t_p->nor_x, t_p->nor_z);
  latitude  = atan2f(t_p->nor_y, sqrtf(t_p->nor_x * t_p->nor_x + t_p->nor_z * t_p->nor_z));
  
  t_p->d.sin_longitude = sinf(longitude);
  t_p->d.cos_longitude = cosf(longitude);
  t_p->d.sin_latitude  = sinf(latitude);
  t_p->d.cos_latitude  = cosf(latitude);
  
  r = visTriangleRadius (t_p);
  
  t_p->d.r2 = r * r;
}


void visEditTriangle (float t_x, float t_y, float longitude, float latitude, float size_factor,
		      int b_id, int t_id, Vis *vis)
{
  float cx, cy, cz;
  float dx1, dy1, dz1;
  float dx2, dy2, dz2;
  float x, y, z;
  
  int i;
  
  Triangle *t_p;
  
  
  t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  visTriangleCenter (t_p, &cx, &cy, &cz);
  
  for (i = 0; i < 3; i++)
    {
      dx1 = t_p->v[i].pos_x - cx;
      dy1 = t_p->v[i].pos_y - cy;
      dz1 = t_p->v[i].pos_z - cz;
      
      visRotate (dx1, dy1, dz1, longitude, latitude, &dx2, &dy2, &dz2);
      
      t_p->v[i].pos_x = size_factor * dx2 + cx;
      t_p->v[i].pos_y = size_factor * dy2 + cy;
      t_p->v[i].pos_z = size_factor * dz2 + cz;
      
      visProject (t_p->v[i].pos_x, t_p->v[i].pos_y, t_p->v[i].pos_z, &x, &y, &z);
      
      x += t_x;
      y += t_y;
      
      visAntiProject (x, y, z, &t_p->v[i].pos_x, &t_p->v[i].pos_y, &t_p->v[i].pos_z);
    }
  visCalculateTriangleData (t_p);
}


void visDeleteTriangle (int b_id, int t_id, Vis *vis)
{
  Boundary *boundary_p = &vis->boundary[ b_id ];
  
  if (t_id == --boundary_p->triangles)
    {
      return;
    }
  memcpy (&boundary_p->triangle[ t_id ],
	  &boundary_p->triangle[ boundary_p->triangles ], sizeof(Triangle));
}


void visInvertTriangleNormal (int b_id, int t_id, Vis *vis)
{
  Triangle *t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  t_p->normal_sign = -t_p->normal_sign;
}


void visChangeTrianglePars (int b_id, int t_id, float dp_avg, float dp_amp, float dp_phs, Vis *vis)
{
  Triangle *t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  t_p->pressure_avg += dp_avg;
  t_p->pressure_avg = fmaxf (0.0, t_p->pressure_avg);
  
  t_p->pressure_amp += dp_amp;
  t_p->pressure_amp = fminf (t_p->pressure_avg, t_p->pressure_amp);
  
  t_p->pressure_phs += dp_phs;
  t_p->pressure_phs = fmaxf (0.0, t_p->pressure_phs);
}


int visSegmentation (Vis *vis)
{
  float seconds = myClock ();
  
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int si, sj, sk;
  int block_id, site_id;
  int sites_b;
  int are_fluid_sites_incrementing;
  int i, l, m;
  
  Site *site_p, *neigh_site_p;
  
  Block *block_p;
  
  Coords *site_coords_p;
  
  
  for (m = 0; m < vis->blocks; m++)
    {
      vis->block[ m ].site = NULL;
      vis->block[ m ].is_void = 1;
    }
  vis->sites = 0;
  vis->boundary_sites = 0;
  vis->stack_sites = 0;
  vis->sites_a = 0;
  vis->iters = 0;
  
  visFromVoxelToSiteCoords (vis->selected_pixel_x, vis->selected_pixel_y, vis->selected_slice,
			    &site_i, &site_j, &site_k, vis);
  
  if (visInterpolatedGray (site_i, site_j, site_k, vis) < vis->selected_gray)
    {
      return !SUCCESS;
    }
  visUpdateTrianglesGrid (vis);
  
  visFromSiteToGridCoords (site_i, site_j, site_k, &si, &sj, &sk, &block_id, &site_id, vis);
  
  (block_p = &vis->block[ block_id ])->site = visStackSitePointer (vis);
  
  for (m = 0; m < vis->sites_in_a_block; m++)
    {
      block_p->site[ m ].data = SOLID_TYPE;
    }
  block_p->site[ site_id ].data = FLUID_TYPE;
  block_p->site[ site_id ].iters = 0;
  block_p->site[ site_id ].label = 0;
  block_p->site[ site_id ].x = si;
  block_p->site[ site_id ].y = sj;
  block_p->site[ site_id ].z = sk;
  block_p->is_void = 0;
  
  vis->sites = 1;
  
  vis->site_coords[ 0 ][ 0 ].x = site_i;
  vis->site_coords[ 0 ][ 0 ].y = site_j;
  vis->site_coords[ 0 ][ 0 ].z = site_k;
  vis->sites_a = 1;

  are_fluid_sites_incrementing = 1;
  
  while (are_fluid_sites_incrementing)
    {
      are_fluid_sites_incrementing = 0;
      ++vis->iters;
      sites_b = 0;
      
      for (i = 0; i < vis->sites_a; i++)
	{
	  site_i = vis->site_coords[ 0 ][ i ].x;
	  site_j = vis->site_coords[ 0 ][ i ].y;
	  site_k = vis->site_coords[ 0 ][ i ].z;
	  
	  site_p = visSitePointer (site_i, site_j, site_k, vis);
	  
	  for (l = 0; l < 14; l++)
	    {
	      neigh_i = site_i + e[ l*3+0 ];
	      neigh_j = site_j + e[ l*3+1 ];
	      neigh_k = site_k + e[ l*3+2 ];
	      
	      if (neigh_i == -1 || neigh_i >= vis->sites_x) continue;
	      if (neigh_j == -1 || neigh_j >= vis->sites_y) continue;
	      if (neigh_k == -1 || neigh_k >= vis->sites_z) continue;
	      
	      visFromSiteToGridCoords (neigh_i, neigh_j, neigh_k, &si, &sj, &sk, &block_id, &site_id, vis);
	      
	      block_p = &vis->block[ block_id ];
	      
	      if (block_p->site != NULL &&
		  (neigh_site_p = &block_p->site[ site_id ])->data != SOLID_TYPE)
		{
		  continue;
		}
	      if (visInterpolatedGray (neigh_i, neigh_j, neigh_k, vis) < vis->selected_gray)
		{
		  continue;
		}
	      if (visIsSegmentIntercepted (site_i, site_j, site_k, l, site_p, vis) == SUCCESS)
		{
		  if ((site_p->data & SITE_TYPE_MASK) == INLET_TYPE ||
		      (site_p->data & SITE_TYPE_MASK) == OUTLET_TYPE)
		    {
		      if (vis->boundary_sites >= COORDS_BUFFERS_SIZE_B)
			{
			  return !SUCCESS;
			}
		      vis->site_coords[ 2 ][ vis->boundary_sites ].x = site_i;
		      vis->site_coords[ 2 ][ vis->boundary_sites ].y = site_j;
		      vis->site_coords[ 2 ][ vis->boundary_sites ].z = site_k;
		      ++vis->boundary_sites;
		    }
		  continue;
		}
	      if (sites_b >= COORDS_BUFFERS_SIZE_A)
		{
		  return !SUCCESS;
		}
	      are_fluid_sites_incrementing = 1;
	      
	      if (block_p->site == NULL)
		{
		  if ((block_p->site = visStackSitePointer (vis)) == NULL)
		    {
		      return !SUCCESS;
		    }
		  for (m = 0; m < vis->sites_in_a_block; m++)
		    {
		      block_p->site[ m ].data = SOLID_TYPE;
		    }
		}
	      neigh_site_p = &block_p->site[ site_id ];
	      neigh_site_p->data = FLUID_TYPE;
	      neigh_site_p->iters = vis->iters;
	      neigh_site_p->label = 0;
	      neigh_site_p->x = si;
	      neigh_site_p->y = sj;
	      neigh_site_p->z = sk;
	      block_p->is_void = 0;
	      
	      ++vis->sites;
	      
	      vis->site_coords[ 1 ][ sites_b ].x = neigh_i;
	      vis->site_coords[ 1 ][ sites_b ].y = neigh_j;
	      vis->site_coords[ 1 ][ sites_b ].z = neigh_k;
	      ++sites_b;
	    }
	}
      site_coords_p = vis->site_coords[ 0 ];
      vis->site_coords[ 0 ] = vis->site_coords[ 1 ];
      vis->site_coords[ 1 ] = site_coords_p;
      vis->sites_a = sites_b;
    }
  vis->segmentation_time = myClock () - seconds;
  
  return SUCCESS;
}


void visCreateSplitTriangle (int b_id, int t_id,
			     float cx, float cy, float cz, float radius, Vis *vis)
{
  Triangle *t1, *t2;
  
  
  t1 = &vis->boundary[ b_id ].triangle[ t_id ];
  t2 = &vis->boundary[ b_id ].triangle[ vis->boundary[b_id].triangles ];
  
  memcpy (t2, t1, sizeof(Triangle));
  
  visRotate (0.F, radius * 2.F, 0.F,
	     t2->d.sin_longitude, t2->d.cos_longitude,
	     t2->d.sin_latitude,  t2->d.cos_latitude,
	     &t2->v[0].pos_x, &t2->v[0].pos_y, &t2->v[0].pos_z);
  
  visRotate (-(radius * sqrtf(3.F)), -radius, 0.F,
	     t2->d.sin_longitude, t2->d.cos_longitude,
	     t2->d.sin_latitude,  t2->d.cos_latitude,
	     &t2->v[1].pos_x, &t2->v[1].pos_y, &t2->v[1].pos_z);
  
  visRotate (+(radius * sqrtf(3.F)), -radius, 0.F,
	     t2->d.sin_longitude, t2->d.cos_longitude,
	     t2->d.sin_latitude,  t2->d.cos_latitude,
	     &t2->v[2].pos_x, &t2->v[2].pos_y, &t2->v[2].pos_z);
  
  for (int l = 0; l < 3; l++)
    {
      t2->v[l].pos_x += cx;
      t2->v[l].pos_y += cy;
      t2->v[l].pos_z += cz;
    }
  t2->d.r2 = radius * radius;
  
  t2->pos_x = cx;
  t2->pos_y = cy;
  t2->pos_z = cz;
  
  ++vis->boundary[ b_id ].triangles;
}


int visOptimiseBoundaries (Vis *vis)
{
  float cx, cy, cz;
  float dx, dy, dz;
  float r;
  
  int is_split[2][(1<<BOUNDARY_ID_BITS)];
  int old_inlet_triangles, old_outlet_triangles;
  int boundary_sites;
  int b_id, t_id;
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int si, sj, sk;
  int block_id, site_id;
  int sites_a, sites_b;
  int are_fluid_sites_incrementing;
  int x, y, z;
  int i, l, n;
  
  Site *site_p, *neigh_site_p;
  
  Block *block_p;
  
  Coords *site_coords_p;
  
  
  old_inlet_triangles  = vis->boundary[ INLET_BOUNDARY ].triangles;
  old_outlet_triangles = vis->boundary[ OUTLET_BOUNDARY ].triangles;
  
  for (n = 0; n < vis->boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      is_split[ 0 ][ n ] = -1;
    }
  for (n = 0; n < vis->boundary[ OUTLET_BOUNDARY ].triangles; n++)
    {
      is_split[ 1 ][ n ] = -1;
    }
  for (int stage_id = 0; stage_id < 2; stage_id++)
    {
      for (n = 0; n < vis->boundary_sites; n++)
	{
	  site_i = vis->site_coords[ 2 ][ n ].x;
	  site_j = vis->site_coords[ 2 ][ n ].y;
	  site_k = vis->site_coords[ 2 ][ n ].z;
	  
	  site_p = visSitePointer (site_i, site_j, site_k, vis);
	  
	  if (site_p->label == 2+stage_id) continue;
	  
	  if ((site_p->data & SITE_TYPE_MASK) == INLET_TYPE)
	    {
	      b_id = INLET_BOUNDARY;
	    }
	  else
	    {
	      b_id = OUTLET_BOUNDARY;
	    }
	  t_id = (site_p->data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
	  
	  memcpy (&vis->site_coords[0][0], &vis->site_coords[2][n], sizeof(Coords));
	  
	  if (stage_id == 1)
	    {
	      if (is_split[ b_id ][ t_id ] <= 0)
		{
		  site_p->label = 2+stage_id;
		  continue;
		}
	      memcpy (&vis->site_coords[3][0], &vis->site_coords[2][n], sizeof(Coords));
	      boundary_sites = 1;
	    }
	  sites_a = 1;
	  are_fluid_sites_incrementing = 1;
	  
	  while (are_fluid_sites_incrementing)
	    {
	      are_fluid_sites_incrementing = 0;
	      sites_b = 0;
	      
	      for (i = 0; i < sites_a; i++)
		{
		  for (l = 0; l < 14; l++)
		    {
		      neigh_i = vis->site_coords[ 0 ][ i ].x + e[ l*3+0 ];
		      neigh_j = vis->site_coords[ 0 ][ i ].y + e[ l*3+1 ];
		      neigh_k = vis->site_coords[ 0 ][ i ].z + e[ l*3+2 ];
		      
		      if (neigh_i == -1 || neigh_i >= vis->sites_x) continue;
		      if (neigh_j == -1 || neigh_j >= vis->sites_y) continue;
		      if (neigh_k == -1 || neigh_k >= vis->sites_z) continue;
		      
		      visFromSiteToGridCoords (neigh_i, neigh_j, neigh_k,
					       &si, &sj, &sk, &block_id, &site_id, vis);
		      
		      if ((block_p = &vis->block[ block_id ])->site == NULL)
			{
			  continue;
			}
		      neigh_site_p = &block_p->site[ site_id ];
		      
		      if (neigh_site_p->label == (2+stage_id) ||
			  (neigh_site_p->data & SITE_TYPE_MASK) != (site_p->data & SITE_TYPE_MASK) ||
			  ((site_p->data & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT) != t_id)
			{
			  continue;
			}
		      if (sites_b >= COORDS_BUFFERS_SIZE_B)
			{
			  return !SUCCESS;
			}
		      are_fluid_sites_incrementing = 1;
		      
		      neigh_site_p->label = 2+stage_id;
		      
		      vis->site_coords[ 1 ][ sites_b ].x = neigh_i;
		      vis->site_coords[ 1 ][ sites_b ].y = neigh_j;
		      vis->site_coords[ 1 ][ sites_b ].z = neigh_k;
		      ++sites_b;
		      

		    }
		}
	      if (are_fluid_sites_incrementing)
		{
		  site_coords_p = vis->site_coords[ 0 ];
		  vis->site_coords[ 0 ] = vis->site_coords[ 1 ];
		  vis->site_coords[ 1 ] = site_coords_p;
		  sites_a = sites_b;
		  
		  if (stage_id == 1)
		    {
		      if (boundary_sites + sites_a >= COORDS_BUFFERS_SIZE_B)
			{
			  return !SUCCESS;
			}
		      memcpy (&vis->site_coords[3][boundary_sites],
			      &vis->site_coords[0][0], sites_a * sizeof(Coords));
		      boundary_sites += sites_a;
		    }
		}
	    }
	  if (stage_id == 0)
	    {
	      ++is_split[ b_id ][ t_id ];
	      continue;
	    }
	  x = 0;
	  y = 0;
	  z = 0;
	  
	  for (i = 0; i < boundary_sites; i++)
	    {
	      x += vis->site_coords[ 3 ][ i ].x;
	      y += vis->site_coords[ 3 ][ i ].y;
	      z += vis->site_coords[ 3 ][ i ].z;
	    }
	  cx = (float)x / boundary_sites;
	  cy = (float)y / boundary_sites;
	  cz = (float)z / boundary_sites;
	  
	  dx = 0.F;
	  dy = 0.F;
	  dz = 0.F;
	  
	  for (i = 0; i < boundary_sites; i++)
	    {
	      dx = fmaxf(dx, vis->site_coords[ 3 ][ i ].x - cx);
	      dy = fmaxf(dy, vis->site_coords[ 3 ][ i ].y - cy);
	      dz = fmaxf(dz, vis->site_coords[ 3 ][ i ].z - cz);
	    }
	  cx -= vis->half_dim_x;
	  cy -= vis->half_dim_y;
	  cz -= vis->half_dim_z;
	  
	  r = 1.1F * sqrtf(dx * dx + dy * dy + dz * dz);
      
	  visCreateSplitTriangle (b_id, t_id, cx, cy, cz, r, vis);
	}
    }
  for (n = 0; n < old_inlet_triangles; n++)
    {
      if (is_split[ 0 ][ n ] == 1) visDeleteTriangle (INLET_BOUNDARY, n, vis);
    }
  for (n = 0; n < old_outlet_triangles; n++)
    {
      if (is_split[ 1 ][ n ] == 1) visDeleteTriangle (OUTLET_BOUNDARY, n, vis);
    }
  return SUCCESS;
}


void visSetBoundaryConfigurations (Vis *vis)
{
  int unknowns;
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int m, n;
  
  unsigned int boundary_config, l;
  
  Site *site_p;
  Site *neigh_site_p;
  
  Block *block_p;
  
  
  n = -1;
  
  for (i = 0; i < vis->blocks_x; i++)
    {
      for (j = 0; j < vis->blocks_y; j++)
	{
	  for (k = 0; k < vis->blocks_z; k++)
	    {
	      if ((block_p = &vis->block[ ++n ])->is_void) continue;
	      
	      m = -1;
	      
	      for (site_i = i * vis->block_size; site_i < i * vis->block_size + vis->block_size; site_i++)
		{
		  for (site_j = j * vis->block_size; site_j < j * vis->block_size + vis->block_size; site_j++)
		    {
		      for (site_k = k * vis->block_size; site_k < k * vis->block_size + vis->block_size; site_k++)
			{
			  if ((site_p = &block_p->site[ ++m ])->data == SOLID_TYPE) continue;
			  
			  boundary_config = 0U;
			  unknowns = 0;
			  
			  for (l = 0U; l < 14U; l++)
			    {
			      neigh_i = site_i + e[ l*3+0 ];
			      neigh_j = site_j + e[ l*3+1 ];
			      neigh_k = site_k + e[ l*3+2 ];
			      
			      if (neigh_i == -1 || neigh_i >= vis->sites_x) continue;
			      if (neigh_j == -1 || neigh_j >= vis->sites_y) continue;
			      if (neigh_k == -1 || neigh_k >= vis->sites_z) continue;
			      
			      neigh_site_p = visSitePointer (neigh_i, neigh_j, neigh_k, vis);
			      
			      if (!(neigh_site_p == NULL || neigh_site_p->data == SOLID_TYPE))
				{
				  boundary_config |= (1U << l);
				}
			      else
				{
				  ++unknowns;
				}
			    }
			  if ((site_p->data & SITE_TYPE_MASK) == FLUID_TYPE)
			    {
			      if (unknowns != 0)
				{
				  site_p->data |= (boundary_config << BOUNDARY_CONFIG_SHIFT);
				}
			    }
			  else
			    {
			      site_p->data |= (boundary_config << BOUNDARY_CONFIG_SHIFT);
			    }
			}
		    }
		}
	    }
	}
    }
}


void visRescaleSystem (Vis *vis)
{
  vis->output_image_pix_x = vis->input_pixels_x * vis->res_factor;
  vis->output_image_pix_y = vis->input_pixels_y * vis->res_factor;
  
  vis->output_slices = (int)(vis->input_slices * (vis->slice_size / vis->pixel_size) *
			     (float)vis->output_image_pix_x / (float)vis->input_pixels_x);
  
  vis->scale_x = vis->res_factor;
  vis->scale_y = vis->res_factor;
  vis->scale_z = vis->res_factor * vis->slice_size / vis->pixel_size;
  
  vis->scale_inv_x = 1.F / vis->scale_x;
  vis->scale_inv_y = 1.F / vis->scale_y;
  vis->scale_inv_z = 1.F / vis->scale_z;
  
  // system parameters setup
  
  vis->block_size = 8;
  vis->shift = 3;
  
  vis->blocks_x = vis->output_image_pix_x >> vis->shift;
  vis->blocks_y = vis->output_image_pix_y >> vis->shift;
  vis->blocks_z = vis->output_slices      >> vis->shift;
  
  if ((vis->blocks_x << vis->shift) < vis->output_image_pix_x) ++vis->blocks_x;
  if ((vis->blocks_y << vis->shift) < vis->output_image_pix_y) ++vis->blocks_y;
  if ((vis->blocks_z << vis->shift) < vis->output_slices     ) ++vis->blocks_z;
  
  vis->sites_x = vis->blocks_x * vis->block_size;
  vis->sites_y = vis->blocks_y * vis->block_size;
  vis->sites_z = vis->blocks_z * vis->block_size;
  
  vis->dim_x = vis->output_image_pix_x;
  vis->dim_y = vis->output_image_pix_y;
  vis->dim_z = vis->output_slices;
  
  vis->half_dim_x = 0.5F * vis->dim_x;
  vis->half_dim_y = 0.5F * vis->dim_y;
  vis->half_dim_z = 0.5F * vis->dim_z;
  
  vis->system_size = fmaxf(vis->dim_x, fmaxf(vis->dim_y, vis->dim_z));
  
  vis->sites_in_a_block = vis->block_size * vis->block_size * vis->block_size;
  
  vis->blocks = vis->blocks_x * vis->blocks_y * vis->blocks_z;
  
  vis->block = (Block *)realloc(vis->block, sizeof(Block) * vis->blocks);
  
  int n = 0;
  
  for (int i = 0; i < vis->sites_x; i+=vis->block_size)
    {
      for (int j = 0; j < vis->sites_y; j+=vis->block_size)
	{
	  for (int k = 0; k < vis->sites_z; k+=vis->block_size)
	    {
	      vis->block[ n ].x = i;
	      vis->block[ n ].y = j;
	      vis->block[ n ].z = k;
	      ++n;
	    }
	}
    }
}


void visRescaleViewpoint (float scale, Vis *vis)
{
  vis->ortho_x *= scale;
  vis->ortho_y *= scale;
  vis->viewpoint_radius *= scale;
  vis->viewport_radius = 0.5F * vis->viewpoint_radius;
  
  visProjection (vis);
}


void visRescaleTriangles (float scale, Vis *vis)
{
  for (int n = 0; n < BOUNDARIES; n++)
    {
      for (int m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  for (int l = 0; l < 3; l++)
	    {
	      vis->boundary[n].triangle[m].v[l].pos_x *= scale;
	      vis->boundary[n].triangle[m].v[l].pos_y *= scale;
	      vis->boundary[n].triangle[m].v[l].pos_z *= scale;
	    }
	  visCalculateTriangleData (&vis->boundary[ n ].triangle[ m ]);
	}
    }
}


int visIsSuperficialSite (int site_i, int site_j, int site_k, Vis *vis)
{
  int neigh_i, neigh_j, neigh_k;
  int l;
  
  Site *site_p;
  
  
  for (l = 0; l < 14; l++)
    {
      neigh_i = site_i + e[ l*3+0 ];
      neigh_j = site_j + e[ l*3+1 ];
      neigh_k = site_k + e[ l*3+2 ];
      
      if (neigh_i == -1 || neigh_i >= vis->sites_x) return SUCCESS;
      if (neigh_j == -1 || neigh_j >= vis->sites_y) return SUCCESS;
      if (neigh_k == -1 || neigh_k >= vis->sites_z) return SUCCESS;
      
      site_p = visSitePointer (neigh_i, neigh_j, neigh_k, vis);
      
      if (site_p == NULL || site_p->data == SOLID_TYPE)
	{
	  return SUCCESS;
	}
    }
  return !SUCCESS;
}


int visEstimateBoundaryNormal (int site_i, int site_j, int site_k, float nor[], Vis *vis)
{
  float org[3];
  float nx, ny, nz;
  float temp;
  
  int neigh_i, neigh_j, neigh_k;
  int si, sj, sk;
  int block_id, site_id;
  int sites_a, sites_b;
  int are_fluid_sites_incrementing;
  int iters;
  int i, l;
  
  Site *site_p;
  
  Block *block_p;
  
  Coords *site_coords_p;
  
  
  org[0] = (float)site_i;
  org[1] = (float)site_j;
  org[2] = (float)site_k;
  
  nor[0] = 0.F;
  nor[1] = 0.F;
  nor[2] = 0.F;
  
  visFromSiteToGridCoords (site_i, site_j, site_k, &si, &sj, &sk, &block_id, &site_id, vis);
  
  block_p = &vis->block[ block_id ];
  
  for (i = 0; i < vis->sites_in_a_block; i++)
    {
      block_p->site[ i ].label = 0;
    }
  block_p->site[ site_id ].label = 1;
  
  vis->site_coords[ 0 ][ 0 ].x = site_i;
  vis->site_coords[ 0 ][ 0 ].y = site_j;
  vis->site_coords[ 0 ][ 0 ].z = site_k;
  sites_a = 1;
  
  are_fluid_sites_incrementing = 1;
  iters = 0;
  
  while (++iters <= ITERS_MAX && are_fluid_sites_incrementing)
    {
      are_fluid_sites_incrementing = 0;
      sites_b = 0;
      
      for (i = 0; i < sites_a; i++)
	{
	  for (l = 0; l < 14; l++)
	    {
	      neigh_i = vis->site_coords[ 0 ][ i ].x + e[ l*3+0 ];
	      neigh_j = vis->site_coords[ 0 ][ i ].y + e[ l*3+1 ];
	      neigh_k = vis->site_coords[ 0 ][ i ].z + e[ l*3+2 ];
	      
	      if (neigh_i == -1 || neigh_i >= vis->sites_x) continue;
	      if (neigh_j == -1 || neigh_j >= vis->sites_y) continue;
	      if (neigh_k == -1 || neigh_k >= vis->sites_z) continue;
	      
	      visFromSiteToGridCoords (neigh_i, neigh_j, neigh_k,
				       &si, &sj, &sk, &block_id, &site_id, vis);
	      
	      if ((block_p = &vis->block[ block_id ])->site == NULL)
		{
		  continue;
		}
	      site_p = &block_p->site[ site_id ];
	      
	      if ((site_p->data & SITE_TYPE_MASK) != FLUID_TYPE ||
		  site_p->label == 1 ||
		  !visIsSuperficialSite (neigh_i, neigh_j, neigh_k, vis))
		{
		  continue;
		}
	      if (sites_b >= COORDS_BUFFERS_SIZE_B)
		{
		  return !SUCCESS;
		}
	      are_fluid_sites_incrementing = 1;
	      
	      site_p->label = 1;
	      
	      vis->site_coords[ 1 ][ sites_b ].x = neigh_i;
	      vis->site_coords[ 1 ][ sites_b ].y = neigh_j;
	      vis->site_coords[ 1 ][ sites_b ].z = neigh_k;
	      ++sites_b;
	      
	      nx = (float)neigh_i - org[0];
	      ny = (float)neigh_j - org[1];
	      nz = (float)neigh_k - org[2];
	      temp = 1.F / sqrtf(nx * nx + ny * ny + nz * nz);
	      nx *= temp;
	      ny *= temp;
	      nz *= temp;
	      
	      nor[0] += nx;
	      nor[1] += ny;
	      nor[2] += nz;
	    }
	}
      site_coords_p = vis->site_coords[ 0 ];
      vis->site_coords[ 0 ] = vis->site_coords[ 1 ];
      vis->site_coords[ 1 ] = site_coords_p;
      sites_a = sites_b;
    }
  temp = 1.F / sqrtf(nor[0] * nor[0] + nor[1] * nor[1] + nor[2] * nor[2]);
  nor[0] *= temp;
  nor[1] *= temp;
  nor[2] *= temp;
  
  return SUCCESS;
}


void visEstimateDiameter (int site_i, int site_j, int site_k, float nor[], float *diameter, Vis *vis)
{
  float segment_length = 0.25F;
  float org[3];
  float x, y, z;
  
  int neigh_i, neigh_j, neigh_k;
  int iters;
  
  Site *neigh_site_p;
  
  
  org[0] = (float)site_i;
  org[1] = (float)site_j;
  org[2] = (float)site_k;
  
  x = org[0] + 2.F * nor[0];
  y = org[1] + 2.F * nor[1];
  z = org[2] + 2.F * nor[2];
  
  neigh_i = (int)x;
  neigh_j = (int)y;
  neigh_k = (int)z;
  
  if (neigh_i < 0 || neigh_i >= vis->sites_x ||
      neigh_j < 0 || neigh_j >= vis->sites_y ||
      neigh_k < 0 || neigh_k >= vis->sites_z)
    {
      nor[0] = -nor[0];
      nor[1] = -nor[1];
      nor[2] = -nor[2];
    }
  else
    {
      neigh_site_p = visSitePointer (neigh_i, neigh_j, neigh_k, vis);
      
      if (neigh_site_p == NULL || neigh_site_p->data == SOLID_TYPE)
	{
	  nor[0] = -nor[0];
	  nor[1] = -nor[1];
	  nor[2] = -nor[2];
	}
    }
  nor[0] *= segment_length;
  nor[1] *= segment_length;
  nor[2] *= segment_length;
  
  x = org[0];
  y = org[1];
  z = org[2];
  
  iters = 0;
  
  neigh_site_p = visSitePointer (site_i, site_j, site_k, vis);
  
  while (neigh_site_p != NULL && neigh_site_p->data != SOLID_TYPE)
    {
      ++iters;
      
      x += nor[0];
      y += nor[1];
      z += nor[2];
      
      neigh_i = (int)x;
      neigh_j = (int)y;
      neigh_k = (int)z;
      
      if (neigh_i < 0 || neigh_i >= vis->sites_x ||
	  neigh_j < 0 || neigh_j >= vis->sites_y ||
	  neigh_k < 0 || neigh_k >= vis->sites_z)
	{
	  neigh_site_p = NULL;
	}
      else
	{
	  neigh_site_p = visSitePointer (neigh_i, neigh_j, neigh_k, vis);
	}
    }
  --iters;
  
  if (iters == 0)
    {
      *diameter = -1.F;
    }
  else
    {
      *diameter = (float)iters * segment_length;
    }
}


int visCreateOptimisedTriangle (unsigned int site_type,
				int site_i, int site_j, int site_k, Vis *vis)
{
  float triangle_factor = 2.F;
  float org[3], nor[3];
  float x, y, z;
  float longitude, latitude;
  float triangle_size;
  float diameter;
  
  int neigh_i, neigh_j, neigh_k;
  int l, m, n;
  int iters, iters_max;
  int longitude_id, latitude_id;
  int t_id;
  
  Site *site_p, *neigh_site_p;
  
  Triangle *t_p;
  
  
  if (vis->boundary[ site_type ].triangles == (int)(1U << BOUNDARY_ID_BITS))
    {
      return -1;
    }
  if (visEstimateBoundaryNormal (site_i, site_j, site_k, nor, vis) == !SUCCESS)
    {
      return -1;
    }
  visEstimateDiameter (site_i, site_j, site_k, nor, &diameter, vis);
  
  if (diameter < 0.F)
    {
      return -1;
    }
  org[0] = (float)site_i;
  org[1] = (float)site_j;
  org[2] = (float)site_k;
  
  // segment center
  org[0] = 0.5F * (org[0] + (org[0] + diameter * nor[0]));
  org[1] = 0.5F * (org[1] + (org[1] + diameter * nor[1]));
  org[2] = 0.5F * (org[2] + (org[2] + diameter * nor[2]));
  
  site_i = (int)org[0];
  site_j = (int)org[1];
  site_k = (int)org[2];
  
  site_p = visSitePointer (site_i, site_j, site_k, vis);
  
  iters_max = 0;
  longitude_id = 0;
  latitude_id = 0;
  
  latitude = 0.F;
  
  for (n = 0; n < 180; n++)
    {
      longitude = 0.F;
      
      for (m = 0; m < 360; m++)
	{
	  visRotate (0.F, 0.F, 1.F,
		     longitude, latitude,
		     &nor[0], &nor[1], &nor[2]);
	  
	  x = org[0];
	  y = org[1];
	  z = org[2];
	  
	  iters = 0;
	  
	  neigh_site_p = site_p;
	  
	  while (neigh_site_p != NULL && neigh_site_p->data != SOLID_TYPE)
	    {
	      ++iters;
	  
	      x += nor[0];
	      y += nor[1];
	      z += nor[2];
	      
	      neigh_i = (int)x;
	      neigh_j = (int)y;
	      neigh_k = (int)z;
	      
	      if (neigh_i < 0 || neigh_i >= vis->sites_x ||
		  neigh_j < 0 || neigh_j >= vis->sites_y ||
		  neigh_k < 0 || neigh_k >= vis->sites_z)
		{
		  neigh_site_p = NULL;
		}
	      else
		{
		  neigh_site_p = visSitePointer (neigh_i, neigh_j, neigh_k, vis);
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
      latitude += 2.F * DEG_TO_RAD;
    }
  if (iters_max == 0)
    {
      return -1;
    }
  longitude = longitude_id * DEG_TO_RAD;
  latitude = latitude_id * 2.F * DEG_TO_RAD;
  
  t_id = vis->boundary[ site_type ].triangles;
  t_p = &vis->boundary[ site_type ].triangle[ t_id ];

  triangle_size = triangle_factor * diameter;
  
  org[0] = org[0] - vis->half_dim_x;
  org[1] = org[1] - vis->half_dim_y;
  org[2] = org[2] - vis->half_dim_z;
  
  visRotate (0.F, triangle_size * 2.F, 0.F,
	     longitude, latitude,
	     &t_p->v[0].pos_x, &t_p->v[0].pos_y, &t_p->v[0].pos_z);
  
  visRotate (-(triangle_size * sqrtf(3.F)), -triangle_size, 0.F,
	     longitude, latitude,
	     &t_p->v[1].pos_x, &t_p->v[1].pos_y, &t_p->v[1].pos_z);
  
  visRotate (+(triangle_size / sqrtf(3.F)), -triangle_size, 0.F,
	     longitude, latitude,
	     &t_p->v[2].pos_x, &t_p->v[2].pos_y, &t_p->v[2].pos_z);
  
  for (l = 0; l < 3; l++)
    {
      t_p->v[l].pos_x += org[0];
      t_p->v[l].pos_y += org[1];
      t_p->v[l].pos_z += org[2];
    }
  t_p->pressure_avg = 80.0F;
  t_p->pressure_amp = 0.0F;
  t_p->pressure_phs = 0.0F;
  
  t_p->normal_sign = 1;
  
  visCalculateTriangleData (t_p);
  
  ++vis->boundary[ site_type ].triangles;
  
  return t_id;
}


ScreenVoxel *visScreenVoxelPointer (int x, int y, Vis *vis)
{
  int voxel_x = (int)(vis->screen_voxels * (float)x / (float)vis->viewport_pixels_x);
  int voxel_y = (int)(vis->screen_voxels * (float)y / (float)vis->viewport_pixels_y);
  
  return &vis->screen_voxel[ voxel_x * vis->screen_voxels + voxel_y ];
}


void visProjectBoundariesToScreenVoxels (Vis *vis)
{
  float x, y, z;
  
  int voxel_i, voxel_j;
  
  ScreenVoxel *voxel_p;
  
  
  for (int n = 0; n < vis->screen_voxels * vis->screen_voxels; n++)
    {
      vis->screen_voxel[ n ].v_z = 0.F;
    }
  for (int n = 0; n < BOUNDARIES; n++)
    {
      for (int m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  Triangle *t_p = &vis->boundary[ n ].triangle[ m ];
	  
	  for (int l = 0; l < 3; l++)
	    {
	      visProject (t_p->v[l].pos_x, t_p->v[l].pos_y, t_p->v[l].pos_z, &x, &y, &z);
	      
	      voxel_i = (int)((vis->screen_voxels / (2.F * screen.max_x)) * (x + screen.max_x));
	      voxel_j = (int)((vis->screen_voxels / (2.F * screen.max_y)) * (y + screen.max_y));
	      
	      if (voxel_i < 0 || voxel_i >= vis->screen_voxels ||
		  voxel_j < 0 || voxel_j >= vis->screen_voxels)
		{
		  continue;
		}
	      voxel_p = &vis->screen_voxel[ voxel_i * vis->screen_voxels + voxel_j ];
	      
	      if (z > voxel_p->v_z)
		{
		  voxel_p->b_id = n;
		  voxel_p->t_id = m;
		  voxel_p->v_id = l;
		  voxel_p->v_z  = z;
		}
	    }
	}
    }
}


void visMoveTriangleVertexWithMouse (int x0, int y0, float z0, Vis *vis)
{
  float x1, y1, z1;
  float x2, y2, z2;
  
  int i, j, k;
  
  
  x1 = screen.max_x * (-1.F + (x0 << 1) / (float)vis->viewport_pixels_x);
  y1 = screen.max_y * (-1.F + (y0 << 1) / (float)vis->viewport_pixels_y);
  z1 = z0;
  
  visAntiProject (x1, y1, z1, &x2, &y2, &z2);
  
  i = vis->mouse.b_id;
  j = vis->mouse.t_id;
  k = vis->mouse.v_id;
  
  vis->boundary[i].triangle[j].v[k].pos_x = x2;
  vis->boundary[i].triangle[j].v[k].pos_y = y2;
  vis->boundary[i].triangle[j].v[k].pos_z = z2;
  
  visCalculateTriangleData (&vis->boundary[i].triangle[j]);
}


void visRotateTriangleWithMouse (float dx1, float dy1,
				 int boundary_index, int triangle_index, Vis *vis)
{
  float cx, cy, cz;
  float nx, ny, nz;
  float dz1;
  float dx2, dy2, dz2;
  float dx3, dy3, dz3;
  float longitude1, latitude1;
  float longitude2, latitude2;
  
  int i;
  
  Triangle *t_p;
  
  
  t_p = &vis->boundary[ boundary_index ].triangle[ triangle_index ];
  
  visTriangleCenter (t_p, &cx, &cy, &cz);
  visTriangleNormal (t_p, &nx, &ny, &nz);
  
  longitude1 = atan2f(nx, nz);
  latitude1  = atan2f(ny, sqrtf(nx * nx + nz * nz));
  
  visRotate (dx1, dy1, 1.F,
	     longitude1, latitude1,
	     &dx2, &dy2, &dz2);
  
  longitude2 = atan2f(dx2, dz2);
  latitude2  = atan2f(dy2, sqrtf(dx2 * dx2 + dz2 * dz2));
  
  for (i = 0; i < 3; i++)
    {
      dx1 = t_p->v[i].pos_x - cx;
      dy1 = t_p->v[i].pos_y - cy;
      dz1 = t_p->v[i].pos_z - cz;
      
      visAntiRotate (dx1, dy1, dz1, longitude1, latitude1, &dx2, &dy2, &dz2);
      
      visRotate (dx2, dy2, dz2, longitude2, latitude2, &dx3, &dy3, &dz3);
      
      t_p->v[i].pos_x = dx3 + cx;
      t_p->v[i].pos_y = dy3 + cy;
      t_p->v[i].pos_z = dz3 + cz;
    }
  visCalculateTriangleData (t_p);
}


void visScaleTriangleWithMouse (float scaling_factor, int b_id, int t_id, Vis *vis)
{
  float cx, cy, cz;
  float dx, dy, dz;
  
  int i;
  
  Triangle *t_p;
  
  
  t_p = &vis->boundary[ b_id ].triangle[ t_id ];
  
  visTriangleCenter (t_p, &cx, &cy, &cz);
  
  for (i = 0; i < 3; i++)
    {
      dx = t_p->v[i].pos_x - cx;
      dy = t_p->v[i].pos_y - cy;
      dz = t_p->v[i].pos_z - cz;
      
      t_p->v[i].pos_x = (1.F + scaling_factor) * dx + cx;
      t_p->v[i].pos_y = (1.F + scaling_factor) * dy + cy;
      t_p->v[i].pos_z = (1.F + scaling_factor) * dz + cz;
    }
  visCalculateTriangleData (t_p);
}


void visRotateViewpointWithMouse (float dx1, float dy1, Vis *vis)
{
  float dx2, dy2, dz2;
  
  
  visRotate (dx1, dy1, 1.F,
	     vis->longitude * DEG_TO_RAD, vis->latitude * DEG_TO_RAD,
	     &dx2, &dy2, &dz2);
  
  vis->longitude = atan2f(dx2, dz2) / DEG_TO_RAD;
  vis->latitude  = atan2f(dy2, sqrtf(dx2 * dx2 + dz2 * dz2)) / DEG_TO_RAD;
  
  visProjection (vis);
}


void visMouseFunction (int button, int state, int x, int y, Vis *vis)
{
  int voxel_i, voxel_j;
  
  
  vis->mouse.state = !ACTIVE;
  
  if (button != GLUT_LEFT_BUTTON) return;
  
  y = vis->viewport_pixels_y - y - 1;
  
  voxel_i = (x * vis->screen_voxels) / vis->viewport_pixels_x;
  voxel_j = (y * vis->screen_voxels) / vis->viewport_pixels_y;
  
  if (voxel_i < 0 || voxel_i >= vis->screen_voxels ||
      voxel_j < 0 || voxel_j >= vis->screen_voxels)
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
  vis->mouse.x = x;
  vis->mouse.y = y;
}


void visMotionFunction (int x, int y, Vis *vis)
{
  if (vis->mouse.state == !ACTIVE) return;
  
  y = vis->viewport_pixels_y - y - 1;
  
  int mouse_dx = x - vis->mouse.x;
  int mouse_dy = y - vis->mouse.y;
  
  if (vis->menu.option == CHANGE_THRESHOLD)
    {
      if (mouse_dy > 0)
	{
	  vis->selected_gray += mouse_dy * mouse_dy;
	}
      else
	{
	  vis->selected_gray -= mouse_dy * mouse_dy;
	}
      vis->selected_gray = fmaxf(vis->gray_min, fminf(vis->gray_max, vis->selected_gray));
    }
  else if (vis->mode == 0)
    {
      if (vis->menu.option == CHANGE_SLICE)
	{
	  vis->selected_slice += mouse_dy;
	  vis->selected_slice = max(0, min(vis->input_slices-1, vis->selected_slice));
	}
    }
  else
    {
      if (vis->menu.option == ZOOM_SCENE)
	{
	  vis->zoom *= 1.F + (float)mouse_dy / vis->viewport_pixels_y;
	  visProjection (vis);
	}
      else if (vis->menu.option == ROTATE_SCENE)
	{
	  visRotateViewpointWithMouse (-(float)mouse_dx / vis->viewport_pixels_x,
				       -(float)mouse_dy / vis->viewport_pixels_y, vis);
	}
      else if (vis->menu.option == ZOOM_BOUNDARY && vis->mouse.b_id >= 0)
	{
	  visScaleTriangleWithMouse ((float)mouse_dy / vis->viewport_pixels_y,
				     vis->mouse.b_id,
				     vis->mouse.t_id,
				     vis);
	}
      else if (vis->menu.option == ROTATE_BOUNDARY && vis->mouse.b_id >= 0)
	{
	  visRotateTriangleWithMouse ((float)mouse_dx / vis->viewport_pixels_x,
				      (float)mouse_dy / vis->viewport_pixels_y,
				      vis->mouse.b_id,
				      vis->mouse.t_id,
				      vis);
	}
    }
  if (vis->mode == 2)
    {
      vis->mode = 1;
    }
  vis->mouse.x = x;
  vis->mouse.y = y;
}


void visPassiveMotionFunction (int x, int y, Vis *vis)
{
  int voxel_i, voxel_j;
  
  ScreenVoxel *voxel_p;
  
  
  y = vis->viewport_pixels_y - y - 1;
  
  vis->mouse.x = x;
  vis->mouse.y = y;
  
  if (vis->mode == 0) return;
  
  voxel_i = x * vis->screen_voxels / vis->viewport_pixels_x;
  voxel_j = y * vis->screen_voxels / vis->viewport_pixels_y;
  
  if (voxel_i < 0 || voxel_i >= vis->screen_voxels ||
      voxel_j < 0 || voxel_j >= vis->screen_voxels)
    {
      vis->mouse.b_id = -1;
      return;
    }
  if (vis->mouse.b_id >= 0) return;
  
  visProjectBoundariesToScreenVoxels (vis);
  
  voxel_p = visScreenVoxelPointer (vis->mouse.x, vis->mouse.y, vis);
  
  if (voxel_p->v_z > EPSILON)
    {
      vis->mouse.t_id = voxel_p->t_id;
      vis->mouse.b_id = voxel_p->b_id;
      vis->mouse.v_id = voxel_p->v_id;
    }
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
  
  visVisualiseString (0.F, 0.F, 0.F, 10, 10,
		      pars_string, GLUT_BITMAP_HELVETICA_12);
}


void visVisualiseSiteData (int site_i, int site_j, int site_k, Vis *vis)
{
  float nor[3];
  float diameter;
  
  char pars_string[256];
  
  
  if (visEstimateBoundaryNormal (site_i, site_j, site_k, nor, vis) == !SUCCESS)
    {
      return;
    }
  visEstimateDiameter (site_i, site_j, site_k, nor, &diameter, vis);
  
  if (diameter >= 0.F)
    {
      diameter *= vis->pixel_size / vis->res_factor;
      sprintf (pars_string, "Lattice coords = (%i, %i, %i)  Diameter = %.3f mm",
	       site_i, site_j, site_k, diameter);
    }
  else
    {
      sprintf (pars_string, "Lattice coords = (%i, %i, %i)  Diameter = NaN",
	       site_i, site_j, site_k);
    }
  visVisualiseString (0.F, 0.F, 0.F, 10, 10,
		      pars_string, GLUT_BITMAP_HELVETICA_12);
}


void visVisualiseVisData (Vis *vis)
{
  char data_string[256];
  
  
  if (vis->mode == 0)
    {
      sprintf (data_string, "Pixel size (mm) = %.3f  Slice thickness (mm) = %.3f  Threshold = %.1f  ",
	       vis->pixel_size, vis->slice_size, vis->selected_gray);
      
      visVisualiseString (0.F, 1.F, 0.F,
			  10, vis->viewport_pixels_y-20,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "Resolution enhancement = %i  Sites = %i  Segmentation time = %.3f",
	       vis->res_factor, vis->sites, vis->segmentation_time);
      
      visVisualiseString (0.F, 1.F, 0.F, 10, vis->viewport_pixels_y-40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
    }
  else
    {
      sprintf (data_string, "Pixel size (mm) = %.3f  Slice thickness (mm) = %.3f  Threshold = %.1f  ",
	       vis->pixel_size, vis->slice_size, vis->selected_gray);
      
      visVisualiseString (0.F, 0.F, 0.F,
			  10, vis->viewport_pixels_y-20,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "Resolution enhancement = %i  Sites = %i  Segmentation time = %.3f",
	       vis->res_factor, vis->sites, vis->segmentation_time);
      
      visVisualiseString (0.F, 0.F, 0.F, 10, vis->viewport_pixels_y-40,
			  data_string, GLUT_BITMAP_HELVETICA_12);
      
      sprintf (data_string, "FPS = %.2f", vis->fps);
  
      visVisualiseString (0.F, 0.F, 0.F,
			  vis->viewport_pixels_x-100, 10,
			  data_string, GLUT_BITMAP_HELVETICA_12);
    }
}


void visVisualiseTriangles (Vis *vis)
{
  float x0, y0, z0;
  float x1, y1, z1;
  float nx, ny, nz;
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
	      glVertex3f (t_p->v[i].pos_x, t_p->v[i].pos_y, t_p->v[i].pos_z);
	    }
	  glEnd ();
	  
	  if (n != INLET_BOUNDARY) continue;
	  
	  glBegin (GL_LINES);
	  
	  x0 = t_p->pos_x;
	  y0 = t_p->pos_y;
	  z0 = t_p->pos_z;
	  
	  nx = t_p->nor_x;
	  ny = t_p->nor_y;
	  nz = t_p->nor_z;

	  length = sqrtf(t_p->d.r2);
	  
	  if (t_p->normal_sign == -1)
	    {
	      nx = -nx;
	      ny = -ny;
	      nz = -nz;
	    }
	  x1 = x0 + nx * length;
	  y1 = y0 + ny * length;
	  z1 = z0 + nz * length;
	  
	  glVertex3f (x0, y0, z0);
	  glVertex3f (x1, y1, z1);
	  
	  glEnd ();
	}
    }
}


void visVisualiseDiscs (Vis *vis)
{
  float v[2*37];
  
  float x1, y1, z1;
  float x2, y2, z2;
  
  int i;
  
  Triangle *t_p;
  
  
  for (i = 0; i <= 35; i++)
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
	  
	  glVertex3f (t_p->pos_x, t_p->pos_y, t_p->pos_z);
	  
	  for (i = 0; i <= 36; i++)
	    {
	      x1 = sqrtf(t_p->d.r2) * v[ 2*i   ];
	      y1 = sqrtf(t_p->d.r2) * v[ 2*i+1 ];
	      
	      visRotate (x1, y1, 0.F,
			 t_p->d.sin_longitude, t_p->d.cos_longitude,
			 t_p->d.sin_latitude,  t_p->d.cos_latitude,
			 &x2, &y2, &z2);
	      
	      x1 = x2 + t_p->pos_x;
	      y1 = y2 + t_p->pos_y;
	      z1 = z2 + t_p->pos_z;
	      
	      glVertex3f (x1, y1, z1);
	    }
	  glEnd ();
	}
    }
}


void visVisualiseActiveBoundaryVoxel (Vis *vis)
{
  float x1, y1;
  float x2, y2;
  
  int voxel_i, voxel_j;
  int matrix_mode;
  
  
  glGetIntegerv (GL_MATRIX_MODE, &matrix_mode);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity();
  
  gluOrtho2D (-screen.max_x, screen.max_x,
	      -screen.max_y, screen.max_y);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix ();
  glLoadIdentity();
  
  glDisable (GL_DEPTH_TEST);
  
  
  voxel_i = (vis->screen_voxels * vis->mouse.x) / vis->viewport_pixels_x;
  voxel_j = (vis->screen_voxels * vis->mouse.y) / vis->viewport_pixels_y;
  
  x1 = screen.max_x * (-1.F + (voxel_i << 1) / (float)vis->screen_voxels);
  y1 = screen.max_y * (-1.F + (voxel_j << 1) / (float)vis->screen_voxels);
  
  x2 = x1 + (screen.max_x * (2.F / vis->screen_voxels));
  y2 = y1 + (screen.max_y * (2.F / vis->screen_voxels));
  
  glColor3f (0.F, 0.F, 1.F);
  glBegin (GL_LINE_LOOP);
  
  glVertex2f (x1, y1);
  glVertex2f (x1, y2);
  glVertex2f (x2, y2);
  glVertex2f (x2, y1);
  
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
  float delta_x, delta_y;
  float gray;
  
  int matrix_mode;
  int i, j, k;
  
  
  glGetIntegerv (GL_MATRIX_MODE, &matrix_mode);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity();
  
  gluOrtho2D (-screen.max_x, screen.max_x,
	      -screen.max_y, screen.max_y);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix ();
  glLoadIdentity();
  
  glDisable (GL_DEPTH_TEST);
  glClear (GL_COLOR_BUFFER_BIT);
  
  delta_x = 2.F * screen.max_x * vis->viewport_pixels_x / (vis->input_pixels_x * (vis->viewport_pixels_x-2));
  delta_y = 2.F * screen.max_y * vis->viewport_pixels_y / (vis->input_pixels_y * (vis->viewport_pixels_y-2));
  
  k = vis->selected_slice;
  
  y = -screen.max_y;
  
  for (j = 0; j < vis->input_pixels_y - 1; j++)
    {
      glBegin (GL_QUAD_STRIP);
      
      x = -screen.max_x;
      
      for (i = 0; i < vis->input_pixels_x; i++)
	{
	  gray = (vis->medical_data[ visVoxelId(i,j,k,vis) ]-vis->gray_min) / (vis->gray_max-vis->gray_min);
	  
	  glColor3f (gray, gray, gray);	  
  	  glVertex2f (x, y);
  	  
	  gray = (vis->medical_data[ visVoxelId(i,j+1,k,vis) ]-vis->gray_min) /(vis->gray_max-vis->gray_min);
	  
	  glColor3f (gray, gray, gray);
  	  glVertex2f (x, y + delta_y);
	  
	  x += delta_x;
  	}
      y += delta_y;
      
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


void visVisualiseFluidSitesFast (Vis *vis)
{
  double seconds = myClock ();
  
  float x1, y1, z1;
  float x2, y2, z2;
  float r, g, b;
  
  int voxel_i, voxel_j;
  int i, j, k;
  int m, n;
  
  ScreenVoxel *screen_voxel_p;
  
  Block *block_p;
  
  
  glShadeModel (GL_FLAT);
  
  glPointSize (2.F);
  glBegin (GL_POINTS);
  
  for (n = 0; n < vis->screen_voxels * vis->screen_voxels; n++)
    {
      vis->screen_voxel[ n ].site_i = -1;
      vis->screen_voxel[ n ].site_z = 0.F;
    }
  for (n = 0; n < vis->blocks; n++)
    {
      block_p = &vis->block[ n ];
      
      if (block_p->site == NULL || block_p->is_void) continue;
      
      for (m = 0; m < vis->sites_in_a_block; m++)
	{
	  if (block_p->site[ m ].data == SOLID_TYPE) continue;
	  
	  x1 = (float)(i = block_p->x + block_p->site[ m ].x) + (0.5F - vis->half_dim_x);
	  y1 = (float)(j = block_p->y + block_p->site[ m ].y) + (0.5F - vis->half_dim_y);
	  z1 = (float)(k = block_p->z + block_p->site[ m ].z) + (0.5F - vis->half_dim_z);
	  
	  visColorPalette (block_p->site[ m ].iters, vis->res_factor, &r, &g, &b);
	  
	  glColor3f (r, g, b);
	  glVertex3f (x1, y1, z1);
	  
	  visProject (x1, y1, z1, &x2, &y2, &z2);
	  
	  voxel_i = (int)((vis->screen_voxels / (2.F * screen.max_x)) * (x2 + screen.max_x));
	  voxel_j = (int)((vis->screen_voxels / (2.F * screen.max_y)) * (y2 + screen.max_y));
	  
	  if (voxel_i < 0 || voxel_i >= vis->screen_voxels ||
	      voxel_j < 0 || voxel_j >= vis->screen_voxels)
	    {
	      continue;
	    }
	  screen_voxel_p = &vis->screen_voxel[ voxel_i * vis->screen_voxels + voxel_j ];
	  
	  if (z2 > screen_voxel_p->site_z)
	    {
	      screen_voxel_p->site_z = z2;
	      screen_voxel_p->site_i = i;
	      screen_voxel_p->site_j = j;
	      screen_voxel_p->site_k = k;
	    }
	}
    }
  glEnd ();
  
  glShadeModel (GL_SMOOTH);
  
  vis->fps = 1.F / (myClock () - seconds);
  
  vis->mode = 2;
}


void visVisualiseFluidSitesSlow (Vis *vis)
{
  double seconds = myClock ();
  
  float x, y, z;
  float r, g, b;
  
  int m, n;
  
  Block *block_p;
  
  
  glShadeModel (GL_FLAT);
  
  for (n = 0; n < vis->blocks; n++)
    {
      block_p = &vis->block[ n ];
      
      if (block_p->site == NULL || block_p->is_void) continue;
      
      for (m = 0; m < vis->sites_in_a_block; m++)
	{
	  if (block_p->site[ m ].data == SOLID_TYPE) continue;
	  
	  x = (float)(block_p->x + block_p->site[ m ].x) - vis->half_dim_x;
	  y = (float)(block_p->y + block_p->site[ m ].y) - vis->half_dim_y;
	  z = (float)(block_p->z + block_p->site[ m ].z) - vis->half_dim_z;
	  
	  visColorPalette (block_p->site[ m ].iters, vis->res_factor, &r, &g, &b);
	  
	  glColor3f (r, g, b);
	  glPushMatrix ();
	  glTranslatef (x, y, z);
	  glCallList (1);
	  glPopMatrix ();
	}
    }
  glShadeModel (GL_SMOOTH);
  
  vis->fps = 1.F / (myClock () - seconds);
}


void visVisualiseSystem (Vis *vis)
{
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  if (vis->mode == 1)
    {
      visVisualiseFluidSitesFast (vis);
    }
  else
    {
      visVisualiseFluidSitesSlow (vis);
    }
  if (vis->mouse.b_id == INLET_BOUNDARY || vis->mouse.b_id == OUTLET_BOUNDARY)
    {
      visVisualiseTrianglePars (vis->mouse.b_id, vis->mouse.t_id, vis);
    }
  else if (vis->mouse.b_id < 0)
    {
      ScreenVoxel *voxel_p = visScreenVoxelPointer (vis->mouse.x, vis->mouse.y, vis);
      
      if (voxel_p->site_i != -1)
	{
	  visVisualiseSiteData (voxel_p->site_i, voxel_p->site_j, voxel_p->site_k, vis);
	}
    }
  visVisualiseVisData (vis);
  visVisualiseTriangles (vis);
  visVisualiseDiscs (vis);
  
  glutSwapBuffers ();
}


void GLUTCALLBACK Visualise (void)
{
  void (*VisualisePointer[2]) (Vis *vis);
  
  VisualisePointer[0] = visVisualiseSelectedSlice;
  VisualisePointer[1] = visVisualiseSystem;
  
  (*VisualisePointer[ min(1, vis.mode) ]) (&vis);
}


int visGetDir (string dir, vector<string> &files)
{
  DIR *dp;
  
  struct dirent *dirp;
  
  
  if ((dp = opendir(dir.c_str())) == NULL)
    {
      cout << "Error(" << errno << ") opening " << dir << endl;
      return errno;
    }
  
  while ((dirp = readdir(dp)) != NULL)
    {
      files.push_back(string(dirp->d_name));
    }
  closedir(dp);
  
  return 0;
}


void visGetFileNames (Vis *vis)
{
  string dir = string(vis->input_path);
  
  vis->file_list = vector<string>();
  
  visGetDir (dir, vis->file_list);
}


void visReadFirstSlice (Vis *vis)
{
  char ppm_type[16];
  
  string file_name;
  
  ifstream input_file;

  
  file_name = string(vis->input_path) + vis->file_list[2];
  
  input_file.open(file_name.c_str());
  
  input_file >> ppm_type;
  input_file >> vis->input_pixels_x;
  input_file >> vis->input_pixels_y;
  input_file >> vis->pixel_depth;
  
  input_file.close();
}


void visReadSlice (int slice_id, Vis *vis)
{
  char ppm_type[16];
  
  string file_name;
  
  ifstream input_file;
  
  
  file_name = string(vis->input_path) + vis->file_list[ slice_id+2 ];
  
  input_file.open(file_name.c_str());
  
  input_file >> ppm_type;
  input_file >> vis->input_pixels_x;
  input_file >> vis->input_pixels_y;
  input_file >> vis->pixel_depth;
  
  for (int j = 0; j < vis->input_pixels_y; j++)
    {
      for (int i = 0; i < vis->input_pixels_x; i++)
	{
	  int temp;
	  
	  input_file >> temp;
	  vis->medical_data[ visVoxelId(i,j,slice_id,vis) ] = temp;
	  
	  vis->gray_min = fminf(vis->gray_min, temp);
	  vis->gray_max = fmaxf(vis->gray_max, temp);
	}
    }
  input_file.close();
}


void visReadConfig (Vis *vis)
{
  visGetFileNames (vis);
  
  vis->input_slices = vis->file_list.size() - 2;
  

  visReadFirstSlice (vis);
  
  vis->output_image_pix_x = vis->input_pixels_x * vis->res_factor;
  vis->output_image_pix_y = vis->input_pixels_y * vis->res_factor;
  
  vis->output_slices = (int)(vis->input_slices * (vis->slice_size / vis->pixel_size) *
			     (float)vis->output_image_pix_x / (float)vis->input_pixels_x);
  
  vis->scale_x = vis->res_factor;
  vis->scale_y = vis->res_factor;
  vis->scale_z = vis->res_factor * vis->slice_size / vis->pixel_size;
  
  vis->scale_inv_x = 1.F / vis->scale_x;
  vis->scale_inv_y = 1.F / vis->scale_y;
  vis->scale_inv_z = 1.F / vis->scale_z;
  
  // system parameters setup
  
  vis->block_size = 8;
  vis->shift = 3;
  
  vis->blocks_x = vis->output_image_pix_x >> vis->shift;
  vis->blocks_y = vis->output_image_pix_y >> vis->shift;
  vis->blocks_z = vis->output_slices      >> vis->shift;
  
  if ((vis->blocks_x << vis->shift) < vis->output_image_pix_x) ++vis->blocks_x;
  if ((vis->blocks_y << vis->shift) < vis->output_image_pix_y) ++vis->blocks_y;
  if ((vis->blocks_z << vis->shift) < vis->output_slices     ) ++vis->blocks_z;
  
  vis->sites_x = vis->blocks_x * vis->block_size;
  vis->sites_y = vis->blocks_y * vis->block_size;
  vis->sites_z = vis->blocks_z * vis->block_size;
  
  vis->dim_x = vis->output_image_pix_x;
  vis->dim_y = vis->output_image_pix_y;
  vis->dim_z = vis->output_slices;
  
  vis->half_dim_x = 0.5F * vis->dim_x;
  vis->half_dim_y = 0.5F * vis->dim_y;
  vis->half_dim_z = 0.5F * vis->dim_z;
  
  vis->system_size = fmaxf(vis->dim_x, fmaxf(vis->dim_y, vis->dim_z));
  
  vis->sites_in_a_block = vis->block_size * vis->block_size * vis->block_size;
  
  vis->blocks = vis->blocks_x * vis->blocks_y * vis->blocks_z;
  
  vis->medical_data = (unsigned short int *)malloc(sizeof(unsigned short int) *
						   vis->input_slices * vis->input_pixels_x * vis->input_pixels_y);
  
  vis->block = (Block *)malloc(sizeof(Block) * vis->blocks);
  
  int n = 0;
  
  for (int i = 0; i < vis->sites_x; i+=vis->block_size)
    {
      for (int j = 0; j < vis->sites_y; j+=vis->block_size)
	{
	  for (int k = 0; k < vis->sites_z; k+=vis->block_size)
	    {
	      vis->block[ n ].site = NULL;
	      vis->block[ n ].x = i;
	      vis->block[ n ].y = j;
	      vis->block[ n ].z = k;
	      ++n;
	    }
	}
    }
  vis->stack_sites_max = vis->sites_in_a_block * 100000;
  
  vis->stack_site = (Site *)malloc(sizeof(Site) * vis->stack_sites_max);
  
  vis->stack_triangles_max = 1000;
  vis->stack_triangle = (StackTriangle *)malloc(sizeof(StackTriangle) * vis->stack_triangles_max);
  
  vis->gray_max = -1.0e9;
  vis->gray_min = +1.0e9;
  
  for (int i = 0; i < vis->input_slices; i++)
    {
      visReadSlice (i, vis);
    }
  vis->site_coords[ 0 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_A);
  vis->site_coords[ 1 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_A);
  vis->site_coords[ 2 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_B);
  vis->site_coords[ 3 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_B);
}


void visWriteCheckpoint (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;
  
  double dummy;
  
  int i, m, n;
  
  unsigned int data;
  
  
  system_config = fopen (vis->checkpoint, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
  
  xdr_float  (&xdr_config, &vis->slice_size);
  xdr_float  (&xdr_config, &vis->pixel_size);
  
  float res_factor = vis->res_factor;
  
  xdr_float  (&xdr_config, &res_factor);
  xdr_int    (&xdr_config, &vis->smoothing_range);
  
  xdr_int    (&xdr_config, &vis->input_pixels_x);
  xdr_int    (&xdr_config, &vis->input_pixels_y);
  xdr_int    (&xdr_config, &vis->input_slices);
  
  xdr_int    (&xdr_config, &vis->selected_pixel_x);
  xdr_int    (&xdr_config, &vis->selected_pixel_y);
  xdr_int   (&xdr_config, &vis->selected_slice);
  xdr_float  (&xdr_config, &vis->selected_gray);
  
  dummy = 1.F;
  
  xdr_double (&xdr_config, &dummy);
  xdr_int    (&xdr_config, &vis->blocks_x);
  xdr_int    (&xdr_config, &vis->blocks_y);
  xdr_int    (&xdr_config, &vis->blocks_z);
  xdr_int    (&xdr_config, &vis->block_size);
  

  for (n = 0; n < vis->input_slices * vis->input_pixels_x * vis->input_pixels_y; n++)
    {
      data = vis->medical_data[ n ];
      
      xdr_u_int (&xdr_config, &data);
    }
  
  for (n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &vis->boundary[ n ].triangles);
      
      for (m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].v[ i ].pos_x);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].v[ i ].pos_y);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].v[ i ].pos_z);
	    }
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].pressure_avg);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].pressure_amp);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].pressure_phs);
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void visReadCheckpoint (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;
  
  double dummy;
  
  float res_factor;
  
  int i, m, n;
  
  unsigned int data;
  
  
  system_config = fopen (vis->checkpoint, "r");
  xdrstdio_create (&xdr_config, system_config, XDR_DECODE);
  
  xdr_float  (&xdr_config, &vis->slice_size);
  xdr_float  (&xdr_config, &vis->pixel_size);
  xdr_float  (&xdr_config, &res_factor);
  
  vis->res_factor = res_factor;
  
  xdr_int    (&xdr_config, &vis->smoothing_range);
  
  xdr_int    (&xdr_config, &vis->input_pixels_x);
  xdr_int    (&xdr_config, &vis->input_pixels_y);
  xdr_int    (&xdr_config, &vis->input_slices);
  
  xdr_int    (&xdr_config, &vis->selected_pixel_x);
  xdr_int    (&xdr_config, &vis->selected_pixel_y);
  xdr_int    (&xdr_config, &vis->selected_slice);
  xdr_float  (&xdr_config, &vis->selected_gray);
  
  xdr_double (&xdr_config, &dummy);
  xdr_int    (&xdr_config, &vis->blocks_x);
  xdr_int    (&xdr_config, &vis->blocks_y);
  xdr_int    (&xdr_config, &vis->blocks_z);
  xdr_int    (&xdr_config, &vis->block_size);
  
  vis->output_image_pix_x = vis->input_pixels_x * vis->res_factor;
  vis->output_image_pix_y = vis->input_pixels_y * vis->res_factor;
  
  vis->output_slices = (int)(vis->input_slices * (vis->slice_size / vis->pixel_size) *
			     (float)vis->output_image_pix_x / (float)vis->input_pixels_x);
  
  vis->scale_x = (float)vis->output_image_pix_x / (float)vis->input_pixels_x;
  vis->scale_y = (float)vis->output_image_pix_y / (float)vis->input_pixels_y;
  vis->scale_z = (float)vis->output_slices / (float)vis->input_slices;
  
  vis->scale_inv_x = 1.F / vis->scale_x;
  vis->scale_inv_y = 1.F / vis->scale_y;
  vis->scale_inv_z = 1.F / vis->scale_z;
  
  vis->sites_x = vis->blocks_x * vis->block_size;
  vis->sites_y = vis->blocks_y * vis->block_size;
  vis->sites_z = vis->blocks_z * vis->block_size;
  
  vis->dim_x = (float)vis->sites_x;
  vis->dim_y = (float)vis->sites_y;
  vis->dim_z = (float)vis->sites_z;
  
  vis->half_dim_x = 0.5F * vis->dim_x;
  vis->half_dim_y = 0.5F * vis->dim_y;
  vis->half_dim_z = 0.5F * vis->dim_z;
  
  vis->system_size = fmaxf(vis->dim_x, fmaxf(vis->dim_y, vis->dim_z));
  
  vis->sites_in_a_block = vis->block_size * vis->block_size * vis->block_size;
  
  vis->blocks = vis->blocks_x * vis->blocks_y * vis->blocks_z;
  
  i = vis->block_size;
  
  vis->shift = 0;

  while (i > 1)
    {
      i >>= 1;

      ++vis->shift;
    }
  
  // initial setup of the medical and system datasets
  
  vis->medical_data = (unsigned short int *)malloc(sizeof(unsigned short int) *
						   vis->input_slices * vis->input_pixels_x * vis->input_pixels_y);
  
  vis->gray_max = -1.0e9;
  vis->gray_min = +1.0e9;
  
  for (n = 0; n < vis->input_slices * vis->input_pixels_x * vis->input_pixels_y; n++)
    {
      xdr_u_int (&xdr_config, &data);
      
      vis->medical_data[ n ] = data;
      
      vis->gray_min = fminf(vis->gray_min, data);
      vis->gray_max = fmaxf(vis->gray_max, data);
    }
  
  vis->block = (Block *)malloc(sizeof(Block) * vis->blocks);
  
  n = 0;
  
  for (int i = 0; i < vis->sites_x; i+=vis->block_size)
    {
      for (int j = 0; j < vis->sites_y; j+=vis->block_size)
	{
	  for (int k = 0; k < vis->sites_z; k+=vis->block_size)
	    {
	      vis->block[ n ].site = NULL;
	      vis->block[ n ].x = i;
	      vis->block[ n ].y = j;
	      vis->block[ n ].z = k;
	      ++n;
	    }
	}
    }
  
  vis->boundary[ INLET_BOUNDARY  ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  vis->boundary[ OUTLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  vis->boundary[ WALL_BOUNDARY   ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  
  for (n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &vis->boundary[ n ].triangles);
      
      for (m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].v[ i ].pos_x);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].v[ i ].pos_y);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].v[ i ].pos_z);
	    }
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].pressure_avg);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].pressure_amp);
	      xdr_float (&xdr_config, &vis->boundary[ n ].triangle[ m ].pressure_phs);
	    }
	}
    }
  xdr_destroy (&xdr_config);
  
  vis->screen_voxel = (ScreenVoxel *)malloc(sizeof(ScreenVoxel) * vis->screen_voxels * vis->screen_voxels);
  
  vis->site_coords[ 0 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_A);
  vis->site_coords[ 1 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_A);
  vis->site_coords[ 2 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_B);
  vis->site_coords[ 3 ] = (Coords *)malloc(sizeof(Coords) * COORDS_BUFFERS_SIZE_B);
  
  visSetSmoothingRange (vis->smoothing_range, vis->res_factor, vis);
  visSegmentation (vis);
}


void visWriteConfig (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;

  double dummy;
  
  int i, j, k;
  int m, n;
  int are_all_solid_sites;
  int flag;
  
  unsigned int site_data;

  Block *block_p;
  
  
  system_config = fopen (vis->output_config, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
  
  dummy = 1.;
  
  xdr_double (&xdr_config, &dummy);
  xdr_int    (&xdr_config, &vis->blocks_x);
  xdr_int    (&xdr_config, &vis->blocks_y);
  xdr_int    (&xdr_config, &vis->blocks_z);
  xdr_int    (&xdr_config, &vis->block_size);
  
  n = -1;
  
  for (i = 0; i < vis->blocks_x; i++)
    {
      for (j = 0; j < vis->blocks_y; j++)
	{
	  for (k = 0; k < vis->blocks_z; k++)
	    {
	      block_p = &vis->block[ ++n ];
	      
	      flag = 0;
	      
	      if (block_p->site == NULL)
		{
		  xdr_int (&xdr_config, &flag);
		  continue;
		}
	      are_all_solid_sites = 1;
	      
	      for (m = 0; m < vis->sites_in_a_block; m++)
		{
		  if (block_p->site[ m ].data != SOLID_TYPE)
		    {
		      are_all_solid_sites = 0;
		      break;
		    }
		}	      
	      if (are_all_solid_sites)
		{
		  xdr_int (&xdr_config, &flag);
		  continue;
		}
	      flag = 1;
	      xdr_int (&xdr_config, &flag);
	      
	      for (m = 0; m < vis->sites_in_a_block; m++)
		{
		  xdr_u_int (&xdr_config, &site_data);
		}
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void visWritePars (Vis *vis)
{
  FILE *pars = fopen (vis->output_pars, "w");
  
  float nx, ny, nz;
  
  int n;
  
  Triangle *t_p;
  
  
  fprintf (pars, "%i\n", vis->boundary[ INLET_BOUNDARY ].triangles);
  
  for (n = 0; n < vis->boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      fprintf (pars, "%f %f %f\n",
	       vis->boundary[ INLET_BOUNDARY ].triangle[n].pressure_avg,
	       vis->boundary[ INLET_BOUNDARY ].triangle[n].pressure_amp,
	       vis->boundary[ INLET_BOUNDARY ].triangle[n].pressure_phs);
    }
  
  fprintf (pars, "%i\n", vis->boundary[ OUTLET_BOUNDARY ].triangles);
  
  for (n = 0; n < vis->boundary[ OUTLET_BOUNDARY ].triangles; n++)
    {
      fprintf (pars, "%f %f %f\n",
	       vis->boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_avg,
	       vis->boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_amp,
	       vis->boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_phs);
    }
  for (n = 0; n < vis->boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      t_p = &vis->boundary[ INLET_BOUNDARY ].triangle[n];
      
      visTriangleNormal (t_p, &nx, &ny, &nz);
      
      if (t_p->normal_sign == 1)
	{
	  nx = -nx;
	  ny = -ny;
	  nz = -nz;
	}
      fprintf (pars, "%f %f %f\n", nx, ny, nz);
    }
  fclose (pars);
}


void visInitBoundaries (Vis *vis)
{
  vis->boundary[ INLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  vis->boundary[ INLET_BOUNDARY ].triangles = 0;
  
  vis->boundary[ OUTLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  vis->boundary[ OUTLET_BOUNDARY ].triangles = 0;
  
  vis->boundary[ WALL_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  vis->boundary[ WALL_BOUNDARY ].triangles = 0;
  
  vis->screen_voxel = (ScreenVoxel *)malloc(sizeof(ScreenVoxel) * vis->screen_voxels * vis->screen_voxels);
}


void visEndBoundaries (Vis *vis)
{
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      free(vis->boundary[ n ].triangle);
      vis->boundary[ n ].triangle = NULL;
      
      vis->boundary[ n ].triangles = 0;
    }
}
 

void visProcessMenuEvents (int option)
{
  int b_id, t_id;
  
  ScreenVoxel *voxel_p;
  
  
  if (option == CHANGE_SLICE)
    {
      vis.menu.option = CHANGE_SLICE;
    }
  else if (option == CHANGE_THRESHOLD)
    {
      vis.menu.option = CHANGE_THRESHOLD;
    }
  else if (option & (SEGMENT_1X|SEGMENT_2X|SEGMENT_3X|SEGMENT_4X|SEGMENT_5X|SEGMENT_6X))
    {
      float res_factor_temp = vis.res_factor;
      
      if      (option == SEGMENT_1X) { vis.res_factor = 1; }
      else if (option == SEGMENT_2X) { vis.res_factor = 2; }
      else if (option == SEGMENT_3X) { vis.res_factor = 3; }
      else if (option == SEGMENT_4X) { vis.res_factor = 4; }
      else if (option == SEGMENT_5X) { vis.res_factor = 5; }
      else if (option == SEGMENT_6X) { vis.res_factor = 6; }
      
      if (vis.mode == 0)
	{
	  vis.selected_pixel_x = (vis.mouse.x * vis.input_pixels_x) / vis.viewport_pixels_x;
	  vis.selected_pixel_y = (vis.mouse.y * vis.input_pixels_y) / vis.viewport_pixels_y;
	  
	  vis.selected_pixel_x = max(0, min(vis.input_pixels_x-1, vis.selected_pixel_x));
	  vis.selected_pixel_y = max(0, min(vis.input_pixels_y-1, vis.selected_pixel_y));
	}
      visSetSmoothingRange (vis.smoothing_range, vis.res_factor, &vis);
      visRescaleSystem (&vis);
      visRescaleTriangles (vis.res_factor / res_factor_temp, &vis);
      
      if (visSegmentation (&vis) == SUCCESS &&
	  visOptimiseBoundaries (&vis) == SUCCESS)
	{
	  visRescaleViewpoint (vis.res_factor / res_factor_temp, &vis);
	}
      else
	{
	  float temp = vis.res_factor;
	  vis.res_factor = res_factor_temp;
	  res_factor_temp = temp;
	  
	  visSetSmoothingRange (vis.smoothing_range, vis.res_factor, &vis);
	  visRescaleSystem (&vis);
	  visRescaleTriangles (vis.res_factor / res_factor_temp, &vis);
	  visSegmentation (&vis);
	  visOptimiseBoundaries (&vis);
	}
    }
  else if (option & (ZOOM_SCENE|ROTATE_SCENE))
    {
      vis.menu.option = option;
      vis.mouse.b_id = -1;
    }
  else if (option & (CREATE_INLET|CREATE_OUTLET|CREATE_WALL))
    {
      voxel_p = visScreenVoxelPointer (vis.mouse.x, vis.mouse.y, &vis);
      
      if (voxel_p->site_i < 0) return;
      
      if      (option == CREATE_INLET)  { b_id = INLET_BOUNDARY; }
      else if (option == CREATE_OUTLET) { b_id = OUTLET_BOUNDARY; }
      else if (option == CREATE_WALL)   { b_id = WALL_BOUNDARY; }
      
      t_id = visCreateOptimisedTriangle (b_id, voxel_p->site_i, voxel_p->site_j, voxel_p->site_k, &vis);
      
      if (t_id == -1) return;
      
      vis.mouse.b_id = b_id;
      vis.mouse.t_id = t_id;
    }
  else if (option & (ZOOM_BOUNDARY|ROTATE_BOUNDARY))
    {
      vis.menu.option = option;
      
      if (vis.mouse.b_id < 0) return;
    }
  else if (option & (REVERSE_INLET_NORMAL|DELETE_BOUNDARY))
    {
      if (vis.mouse.b_id < 0) return;
      
      if (option == REVERSE_INLET_NORMAL && vis.mouse.b_id == INLET_BOUNDARY)
	{
	  visInvertTriangleNormal (vis.mouse.b_id, vis.mouse.t_id, &vis);
	}
      else if (option == DELETE_BOUNDARY)
	{
	  visDeleteTriangle (vis.mouse.b_id, vis.mouse.t_id, &vis);
	  
	  vis.mouse.b_id = -1;
	}
    }
  else if (option == CHANGE_VIS_MODE)
    {
      vis.menu.option = NULL_MENU_OPTION;
      
      if (vis.mode == 0)
	{
	  vis.mode = 1;
	  glutChangeToMenuEntry (8,  "Zoom scene", ZOOM_SCENE);
	  glutChangeToMenuEntry (9,  "Rotate scene", ROTATE_SCENE);
	  glutChangeToMenuEntry (10, "Create inlet", CREATE_INLET);
	  glutChangeToMenuEntry (11, "Create outlet", CREATE_OUTLET);
	  glutAddMenuEntry ("Create wall", CREATE_WALL);
	  glutAddMenuEntry ("Scale boundary", ZOOM_BOUNDARY);
	  glutAddMenuEntry ("Rotate boundary", ROTATE_BOUNDARY);
	  glutAddMenuEntry ("Reverse inlet normal", REVERSE_INLET_NORMAL);
	  glutAddMenuEntry ("Delete boundary", DELETE_BOUNDARY);
	  glutAddMenuEntry ("2D rendering", CHANGE_VIS_MODE);
	  glutAddMenuEntry ("Save data", SAVE_DATA);
	  glutAddMenuEntry ("Quit", QUIT);
	}
      else
	{
	  vis.mode = 0;
	  vis.mouse.state = !ACTIVE;
	  vis.mouse.b_id = -1;
	  
	  glutChangeToMenuEntry (8,  "Change slice", CHANGE_SLICE);
	  glutChangeToMenuEntry (9,  "3D rendering", CHANGE_VIS_MODE);
	  glutChangeToMenuEntry (10, "Save data", SAVE_DATA);
	  glutChangeToMenuEntry (11, "Quit", QUIT);
	  
	  for (int i = glutGet(GLUT_MENU_NUM_ITEMS); i >= 12; i--)
	    {
	      glutRemoveMenuItem (i);
	    }
	}
    }
  else if (option == SAVE_DATA)
    {
      printf("Opening ppm file: ./image.ppm\n");
      visSaveWindowImage ("./image.ppm");
      
      visSetBoundaryConfigurations (&vis);
      
      printf("Opening output pars file: %s\n", vis.output_pars);
      visWritePars (&vis);
      
      printf("Opening output config file: %s\n", vis.output_config);
      visWriteConfig (&vis);
      
      printf("Opening checkpoint file: %s\n", vis.checkpoint);
      visWriteCheckpoint (&vis);
    }
  else if (option == QUIT)
    {
      visEndBoundaries (&vis);
      
      visEnd (&vis);
      
      exit(0);
    }
}


void visCreateMenu (Vis *vis)
{
  vis->menu.option = NULL_MENU_OPTION;
  vis->mouse.state = !ACTIVE;
  vis->mouse.b_id = -1;
  
  vis->menu.id = glutCreateMenu (visProcessMenuEvents);
  
  glutAddMenuEntry ("Segmentation 1X", SEGMENT_1X);
  glutAddMenuEntry ("Segmentation 2X", SEGMENT_2X);
  glutAddMenuEntry ("Segmentation 3X", SEGMENT_3X);
  glutAddMenuEntry ("Segmentation 4X", SEGMENT_4X);
  glutAddMenuEntry ("Segmentation 5X", SEGMENT_5X);
  glutAddMenuEntry ("Segmentation 6X", SEGMENT_6X);
  glutAddMenuEntry ("Change threshold", CHANGE_THRESHOLD);
  glutAddMenuEntry ("Change slice", CHANGE_SLICE);
  glutAddMenuEntry ("3D rendering", CHANGE_VIS_MODE);
  glutAddMenuEntry ("Save data", SAVE_DATA);
  glutAddMenuEntry ("Quit", QUIT);
  glutAttachMenu (GLUT_RIGHT_BUTTON);
}


void visUsage (char *progname)
{
  printf ("Usage: %s input path, output config, pars and checkpoint names\n", progname);
  printf ("slice and pixel size (mm)\n");
  printf ("    or\n");
  printf ("checkpoint file name, output config and pars file names\n");
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
      
      vis->selected_slice = 0;
      vis->selected_gray = 3500;
      
      vis->screen_voxels = 100;
      vis->mode = 0;
      
      visSetSmoothingRange (1, vis->res_factor, vis);
      
      visReadConfig (vis);
      
      visInitBoundaries (vis);
    }
  else
    {
      vis->checkpoint    = argv[1];
      vis->output_config = argv[2];
      vis->output_pars   = argv[3];
      
      vis->screen_voxels = 100;
      vis->mode = 1;
      
      visReadCheckpoint (vis);
    }
  
  vis->viewport_pixels_x = 512;
  vis->viewport_pixels_y = 512;
  
  vis->background_r = 1.F;
  vis->background_g = 1.F;
  vis->background_b = 1.F;
  
  vis->longitude = 90.F;
  vis->latitude = 0.F;
  
  vis->zoom = 1.0F;
  
  vis->scene_center_x = 0.F;
  vis->scene_center_y = 0.F;
  vis->scene_center_z = 0.F;
  
  vis->ortho_x = 0.5F * vis->system_size;
  vis->ortho_y = 0.5F * vis->system_size;
  vis->viewpoint_radius = 2.F * vis->system_size;
  vis->viewport_radius = 0.5F * vis->viewpoint_radius;
  
  vis->sites = 0;
  vis->segmentation_time = 0.F;
  
  vis->point_size = 2;
  vis->visualise_boundaries = 1;
  
  glutInit (&argc, argv);
  
  visOpenWindow (vis->viewport_pixels_x, vis->viewport_pixels_y);
  
  visProjection (vis);
  
  visCreateCubeDisplayList ();
  
  visCreateMenu (vis);
}


void visEnd (Vis *vis)
{
  visDeleteCubeDisplayList ();
  
  free(vis->screen_voxel);
  
  for (int i = 0; i < COORDS_BUFFERS; i++)
    {
      free(vis->site_coords[ i ]);
    }
  
  free(vis->stack_site);
  
  free(vis->block);
  
  free(vis->medical_data);
}


void visMoveSceneCenter (float t_x, float t_y, Vis *vis)
{
  float x, y, z;
  
  
  visProject (vis->scene_center_x, vis->scene_center_y, vis->scene_center_z, &x, &y, &z);
      
  x += t_x;
  y += t_y;
  
  visAntiProject (x, y, z, &vis->scene_center_x, &vis->scene_center_y, &vis->scene_center_z);
}


void GLUTCALLBACK KeybordFunction (unsigned char key, int x, int y)
{
  if (key == 'c')
    {
      vis.ortho_x = 0.5F * vis.system_size;
      vis.ortho_y = 0.5F * vis.system_size;
      
      vis.longitude = 45.F;
      vis.latitude = 45.F;
      vis.viewpoint_radius = 2.F * vis.system_size;
      vis.viewport_radius = 0.5F * vis.viewpoint_radius;
      
      vis.zoom = 1.0F;
      
      vis.scene_center_x = 0.F;
      vis.scene_center_y = 0.F;
      vis.scene_center_z = 0.F;
      
      visProjection (&vis);
    }
  else if (key == 'j')
    {
      if (vis.mouse.b_id == INLET_BOUNDARY || vis.mouse.b_id == OUTLET_BOUNDARY)
	{
	  visChangeTrianglePars (vis.mouse.b_id, vis.mouse.t_id, 1.0, 0., 0., &vis);
	}
    }
  else if (key == 'J')
    {
      if (vis.mouse.b_id == INLET_BOUNDARY || vis.mouse.b_id == OUTLET_BOUNDARY)
	{
	  visChangeTrianglePars (vis.mouse.b_id, vis.mouse.t_id, -1.0, 0., 0., &vis);
	}
    }
  else if (key == 'k')
    {
      if (vis.mouse.b_id == INLET_BOUNDARY || vis.mouse.b_id == OUTLET_BOUNDARY)
	{
	  visChangeTrianglePars (vis.mouse.b_id, vis.mouse.t_id, 0., 0.1, 0., &vis);
	}
    }
  else if (key == 'K')
    {
      if (vis.mouse.b_id == INLET_BOUNDARY || vis.mouse.b_id == OUTLET_BOUNDARY)
	{
	  visChangeTrianglePars (vis.mouse.b_id, vis.mouse.t_id, 0., -0.1, 0., &vis);
	}
    }
  else if (key == 'l')
    {
      if (vis.mouse.b_id == INLET_BOUNDARY || vis.mouse.b_id == OUTLET_BOUNDARY)
	{
	  visChangeTrianglePars (vis.mouse.b_id, vis.mouse.t_id, 0., 0., 1., &vis);
	}
    }
  else if (key == 'L')
    {
      if (vis.mouse.b_id == INLET_BOUNDARY || vis.mouse.b_id == OUTLET_BOUNDARY)
	{
	  visChangeTrianglePars (vis.mouse.b_id, vis.mouse.t_id, 0., 0., -1., &vis);
	}
    }
}


void GLUTCALLBACK MouseFunction (int button, int state, int x, int y)
{
  visMouseFunction (button, state, x, y, &vis);
}


void GLUTCALLBACK MotionFunction (int x, int y)
{
  visMotionFunction (x, y, &vis);
}


void GLUTCALLBACK PassiveMotionFunction (int x, int y)
{
  visPassiveMotionFunction (x, y, &vis);
}


void GLUTCALLBACK Reshape (GLsizei w, GLsizei h)
{
  // the window is reshaped if necessary
  
  vis.ortho_x *= (float)w / (float)vis.viewport_pixels_x;
  vis.ortho_y *= (float)h / (float)vis.viewport_pixels_y;
  
  vis.viewport_pixels_x = w;
  vis.viewport_pixels_y = h;
  
  glViewport(0, 0, w, h);
  
  visProjection (&vis);
}


int main (int argc, char *argv[])
{
  visInit (argc, argv, &vis);
  
  
  glutReshapeFunc (Reshape);
  glutIdleFunc (Visualise);
  glutDisplayFunc (Visualise);
  glutKeyboardFunc (KeybordFunction);
  glutMouseFunc (MouseFunction);
  glutMotionFunc (MotionFunction);
  glutPassiveMotionFunc (PassiveMotionFunction);
  glutMainLoop ();

  return(0);
}
