#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include <rpc/xdr.h>

#include <iostream>
#include <fstream>

using namespace std;


#ifndef GLUTCALLBACK
#define GLUTCALLBACK
#endif

#define BUFFERS_SIZE   1000000


float EPSILON = 1.0e-30F;
float DEG_TO_RAD = 0.01745329F;

unsigned int SOLID_TYPE  = 0U;
unsigned int FLUID_TYPE  = 1U;
unsigned int INLET_TYPE  = 2U;
unsigned int OUTLET_TYPE = 3U;
unsigned int NULL_TYPE   = 4U;

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

int e_x[] = { 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
int e_y[] = { 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1};
int e_z[] = { 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1};

int inv_dir[] = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12};



char *input_path;
char *output_config_name;
char *output_pars_name;
char *checkpoint_name;

char ppm_type[20];


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

Screen screen;


struct Viewpoint
{
  float pos_x;
  float pos_y;
  float pos_z;
  
  float sin_1;
  float cos_1;
  
  float sin_2;
  float cos_2;
};

Viewpoint viewpoint;


struct LinkedTriangle
{
  short int boundary_id;
  short int triangle_id;
  
  LinkedTriangle *next;
};


struct DataBlock
{
  unsigned int *site_data;
  
  short int *site_label;
  
  unsigned short int *site_iters;
  
  LinkedTriangle *linked_triangle;
  
  int is_void;
};


struct Vertex
{
  float pos_x, pos_y, pos_z;
};


struct Triangle
{
  Vertex v[3];
  
  unsigned int boundary_id, boundary_dir, boundary_config;
  
  float pressure_avg, pressure_amp, pressure_phs;
};


struct Boundary
{
  Triangle *triangle;
  
  int triangles;
};


struct ScreenVoxel
{
  float site_z;
  float vertex_z;
  
  short int site_i, site_j, site_k;
  
  short int triangle_id;
  
  char boundary_id;
  char vertex_id;
};


struct ScreenVoxelCoords
{
  int i, j;
};


struct LastTriangle
{
  int triangle_id;
  int boundary_id;
};


struct SiteLocation
{
  short int i, j, k;
};

struct BlockLocation
{
  short int i, j, k;
};

struct SuperficialSite
{
  short int i, j, k;
  short int iters;
};


struct MyGL
{
  float scale_x, scale_y, scale_z;
  float scale_inv_x, scale_inv_y, scale_inv_z;
  
  int input_image_pix_x, input_image_pix_y;
  int output_image_pix_x, output_image_pix_y;
  int input_slices, output_slices;
  int sites_x, sites_y, sites_z;
  int blocks_x, blocks_y, blocks_z;
  int blocks;
  int block_size;
  int shift;
  int sites_in_a_block;
  int superficial_sites, superficial_sites_max;
  int stored_blocks, stored_blocks_max;
  
  int selected_pixel_x, selected_pixel_y, selected_slice;
  
  unsigned short int selected_gray;
  
  float system_size;
  float dim_x, dim_y, dim_z;
  float half_dim_x, half_dim_y, half_dim_z;
  float lattice_to_system;
  
  DataBlock *data_block;
  
  unsigned short int *medical_data;
  void *slice_row_data;
  
  SuperficialSite *superficial_site;
  
  Boundary boundary[4];
  
  ScreenVoxel *screen_to_boundaries_map;
  
  ScreenVoxel *screen_voxel_pointed;
  
  ScreenVoxelCoords screen_voxel_coords;
  
  int *stored_block;
};


MyGL mygl;


LastTriangle last_triangle;


SiteLocation *site_location_a;
SiteLocation *site_location_b;


int passive_mouse_pixel_i = -1;
int passive_mouse_pixel_j = -1;


int viewport_pixels_x = 1024, viewport_pixels_y = 1024;

float background_r = 1.F, background_g = 1.F, background_b = 1.F;

float ortho_x, ortho_y;

float longitude = 0.F, latitude = 0.F;
float viewpoint_radius;

float zoom = 1.5F;

float scene_center_x = 0.F;
float scene_center_y = 0.F;
float scene_center_z = 0.F;

int screen_voxels = 100;

int point_size = 2;
int draw_triangles = 1;

int display_id;

float screen_voxels_screen_max_inv_x, screen_voxels_screen_max_inv_y;


float slice_size, pixel_size;
float res_factor, res_factor_par1, res_factor_par2;

int voxel_di, voxel_dj, voxel_dk;


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


float myClock ()
{
  return (float)clock () * (1. / (float)CLOCKS_PER_SEC);
}


void myglDisplayString (int x, int y, char *string, void *font)
{
  char *allowed_char[] = {",",".","0","1","2","3","4","5","6","7","8","9"};
  int length = (int)strlen(string);
  
  
  glColor3f (0.F, 0.F, 0.F);
  glRasterPos2f (x, y);
  
  
  for (int i = 0; i < length; i++)
    {
      for (int j = 0; j < 12; j++)//sizeof(allowed_char) / sizeof(char)
  	{
  	  if (!strncmp (&string[i], allowed_char[j],1))
  	    {
  	      glutBitmapCharacter (font, string[i]);
  	    }
  	}
    }
  glEnd ();
}


void myglTriangleNormal (float  x0, float  y0, float  z0,
			 float  x1, float  y1, float  z1,
			 float  x2, float  y2, float  z2,
			 float *nx, float *ny, float *nz)
{
  float dx1, dy1, dz1;
  float dx2, dy2, dz2;
  
  float temp;
  
  
  dx1 = x1 - x0;
  dy1 = y1 - y0;
  dz1 = z1 - z0;
  dx2 = x2 - x0;
  dy2 = y2 - y0;
  dz2 = z2 - z0;
  
  *nx = dy2 * dz1 - dz2 * dy1;
  *ny = dz2 * dx1 - dx2 * dz1;
  *nz = dx2 * dy1 - dy2 * dx1;
  
  temp = 1.F / fmaxf(1.0e-30F, sqrtf(*nx * *nx + *ny * *ny + *nz * *nz));
  
  *nx *= temp;
  *ny *= temp;
  *nz *= temp;
}


float myglTriangleArea (float x1, float y1, float z1,
			float x2, float y2, float z2,
			float x3, float y3, float z3)
{
  float vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3;
  float a, b, c, p;
  
  
  vx1 = x2 - x1;
  vy1 = y2 - y1;
  vz1 = z2 - z1;
  vx2 = x3 - x1;
  vy2 = y3 - y1;
  vz2 = z3 - z1;
  vx3 = x3 - x2;
  vy3 = y3 - y2;
  vz3 = z3 - z2;
  
  a = sqrtf(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);
  b = sqrtf(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);
  c = sqrtf(vx3 * vx3 + vy3 * vy3 + vz3 * vz3);
  
  p = 0.5F * (a + b + c);
  
  return sqrtf(p * (p - a) * (p - b) * (p - c));
}


int myglRayTriangleIntersection (float  px, float  py, float  pz,
				 float  nx, float  ny, float  nz,
				 float vx1, float vy1, float vz1,
				 float vx2, float vy2, float vz2,
				 float vx3, float vy3, float vz3,
				 float  *t, float  *v, float  *w)
{
  float ex1, ey1, ez1, ex2, ey2, ez2;
  float det;
  float x1, y1, z1, x2, y2, z2, x3, y3, z3;
  
  
  ex1 = vx2 - vx1;
  ey1 = vy2 - vy1;
  ez1 = vz2 - vz1;
  
  ex2 = vx3 - vx1;
  ey2 = vy3 - vy1;
  ez2 = vz3 - vz1;
  
  x1 = ny * ez2 - nz * ey2;
  y1 = nz * ex2 - nx * ez2;
  z1 = nx * ey2 - ny * ex2;
  
  det = ex1 * x1 + ey1 * y1 + ez1 * z1;
  
  if (det > -1.e-3F && det < 1.e-30F)
    {
      return 0;
    }
  det = 1.F / det;
  
  x2 = px - vx1;
  y2 = py - vy1;
  z2 = pz - vz1;
  
  *v = (x2 * x1 + y2 * y1 + z2 * z1) * det;
  
  if (*v < 0.F || *v > 1.F)
    {
      return 0;
    }
  
  x3 = y2 * ez1 - z2 * ey1;
  y3 = z2 * ex1 - x2 * ez1;
  z3 = x2 * ey1 - y2 * ex1;
  
  *w = (nx * x3 + ny * y3 + nz * z3) * det;
  
  if (*w < 0.F || *v + *w > 1.F)
    {
      return 0;
    }

  if ((*t = (ex2 * x3 + ey2 * y3 + ez2 * z3) * det) < EPSILON)
    {
      return 0;
    }
  
  return 1;
}


void myglOpenWindow (int pixels_x, int pixels_y)
{
  screen.pix_x = pixels_x;
  screen.pix_y = pixels_y;
  
  glutInitDisplayMode (GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutInitWindowPosition (0, 0);
  glutInitWindowSize (pixels_x, pixels_y);
  
  glutCreateWindow (" ");
  
  glEnable (GL_DEPTH_TEST);
  glDisable (GL_BLEND);
  glShadeModel (GL_FLAT);
  glDisable(GL_DITHER);
  
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}


void myglProjection (float ortho_r, float ortho_g, float ortho_b,
		     float ortho_x, float ortho_y,
		     int pixels_x, int pixels_y,
		     float ctr_x, float ctr_y, float ctr_z,
		     float rad,
		     float longitude, float latitude,
		     float dist,
		     float zoom)
{
  float temp;
  
  
  screen.col_r = ortho_r;
  screen.col_g = ortho_g;
  screen.col_b = ortho_b;
  
  screen.max_x = ortho_x / zoom;
  screen.max_y = ortho_y / zoom;
  
  screen.pix_x = pixels_x;
  screen.pix_y = pixels_y;
  
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
  
  glLoadIdentity();
  glOrtho (-screen.max_x, screen.max_x,
	   -screen.max_y, screen.max_y,
	   -1.F, 1.F);
  
  glClearColor (ortho_r, ortho_g, ortho_b, 0.F);
  
  screen_voxels_screen_max_inv_x = screen_voxels / (2.F * screen.max_x);
  screen_voxels_screen_max_inv_y = screen_voxels / (2.F * screen.max_y);
}


void myglTransformVertex (float  px1, float  py1, float  pz1,
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
  
  *px2 = x2;
  *py2 = y2;

  if (z2 < 0.F)
    {
      z2 = -z2;
      *pz2 = 1.F / (1.F + z2);
    }
  else
    {
      *pz2 = 1.F;
    }
}


void myglAntiTransformVertex (float  px1, float  py1, float  pz1,
			      float *px2, float *py2, float *pz2)
{
  float x1, y1, z1, x2, y2, z2;
  float temp;
  
  
  x1 = px1;
  y1 = py1;
  z1 = 1.F - 1.F / pz1;

  temp = viewpoint.cos_2 * z1 - viewpoint.sin_2 * y1;

  x2 = viewpoint.sin_1 * temp + viewpoint.cos_1 * x1;
  y2 = viewpoint.sin_2 * z1   + viewpoint.cos_2 * y1;
  z2 = viewpoint.cos_1 * temp - viewpoint.sin_1 * x1;

  *px2 = x2 + viewpoint.pos_x;
  *py2 = y2 + viewpoint.pos_y;
  *pz2 = z2 + viewpoint.pos_z;
}


void SaveWindowImage (char *file_name)
{
  FILE *ppm_image_file_ptr = fopen (file_name, "wb");
  
  int pix_x, pix_y;
  int i, j;
  
  unsigned char *data = NULL;  
  unsigned char *data_p = NULL;
  unsigned char *buffer;
  
  
  glReadBuffer (GL_FRONT);
  
  pix_x = screen.pix_x;
  pix_y = screen.pix_y;
  
  data = (unsigned char *)malloc(sizeof(unsigned char) * pix_x * pix_y * 3);
  
  buffer = (unsigned char *)malloc(sizeof(unsigned char) * pix_x * 3);
  
  data_p = data;

  for (j = 0; j < pix_y; j++)
    {
      glReadPixels (0, j, pix_x, 1, GL_RGB, GL_UNSIGNED_BYTE, buffer);
      
      for (i = 0; i < pix_x; i++)
	{
	  *data_p = buffer[ i * 3     ]; data_p++;
	  *data_p = buffer[ i * 3 + 1 ]; data_p++;
	  *data_p = buffer[ i * 3 + 2 ]; data_p++;
	}
    }
  
  free((unsigned char *)buffer);
  
  fprintf (ppm_image_file_ptr, "P6\n%i %i\n255\n", pix_x, pix_y);
  
  for (j = pix_y - 1; j >= 0; j--)
    {
      fwrite (data + j * pix_x * 3, 1, pix_x * 3, ppm_image_file_ptr);
    }
  
  free((unsigned char *)data);
  
  fclose (ppm_image_file_ptr);
}


void Projection (void)
{
    myglProjection (background_r, background_g, background_b,
		    ortho_x, ortho_y,
		    viewport_pixels_x, viewport_pixels_y,
		    scene_center_x, scene_center_y, scene_center_z,
		    viewpoint_radius,    // distance of the viewpoint from the scene center
		    longitude, latitude,    // polar and vertical angles defined by the viewpoint position
		    0.5F * viewpoint_radius,    // distance of the screen from the viewpoint position
		    zoom);
}


float RedComponent (int iters)
{
  iters = iters%400;
  
  if (iters > 200) iters = 400 - iters;
  
  return 0.005F * (float)iters;
}


int myglVoxelId (int i, int j, int k)
{
  return (i * mygl.input_image_pix_y + j) * mygl.input_slices + k;
}


void myglFromVoxelToSiteCoords (int pixel_i, int pixel_j, int slice_id,
				int *site_i, int *site_j, int *site_k)
{
  *site_i = (int)(pixel_i  * mygl.scale_x);
  *site_j = (int)(pixel_j  * mygl.scale_y);
  *site_k = (int)(slice_id * mygl.scale_z);
}


void myglFromSiteToVoxelCoords (int site_i, int site_j, int site_k,
				int *pixel_i, int *pixel_j, int *slice_id)
{
  *pixel_i  = (int)(site_i * mygl.scale_inv_x);
  *pixel_j  = (int)(site_j * mygl.scale_inv_y);
  *slice_id = (int)(site_k * mygl.scale_inv_z);
}


void myglFromSiteToGridCoords (int site_i, int site_j, int site_k, int *block_id, int *site_id)
{
  int i, j, k;
  int ii, jj, kk;
  
  
  i = site_i >> mygl.shift;
  j = site_j >> mygl.shift;
  k = site_k >> mygl.shift;
  *block_id = (i * mygl.blocks_y + j) * mygl.blocks_z + k;
  
  ii = site_i - (i << mygl.shift);
  jj = site_j - (j << mygl.shift);
  kk = site_k - (k << mygl.shift);
  *site_id = (((ii << mygl.shift) + jj) << mygl.shift) + kk;
}

/*
unsigned short int myglInterpolatedGray (int site_i, int site_j, int site_k)
{
  float gray[2][2][2];
  float x, y, z;
  float interpolated_gray;
  
  int voxel_i[2], voxel_j[2], voxel_k[2];
  int i, j, k;
  
  
  x = (float)site_i * mygl.scale_inv_x;
  y = (float)site_j * mygl.scale_inv_y;
  z = (float)site_k * mygl.scale_inv_z;
  
  voxel_i[0] = (int)x;
  voxel_j[0] = (int)y;
  voxel_k[0] = (int)z;
  
  x -= (float)voxel_i[0];
  y -= (float)voxel_j[0];
  z -= (float)voxel_k[0];
  
  voxel_i[1] = min(voxel_i[0] + 1, mygl.input_image_pix_x - 1);
  voxel_j[1] = min(voxel_j[0] + 1, mygl.input_image_pix_y - 1);
  voxel_k[1] = min(voxel_k[0] + 1, mygl.input_slices - 1);
  
  for (i = 0; i < 2; i++)
    {
      for (j = 0; j < 2; j++)
	{
	  for (k = 0; k < 2; k++)
	    {
	      gray[i][j][k] = (float)mygl.medical_data[ myglVoxelId(voxel_i[i],voxel_j[j],voxel_k[k]) ];
	    }
	}
    }
  float gray_00z = (1.F - z) * gray[0][0][0] + z * gray[0][0][1];
  float gray_01z = (1.F - z) * gray[0][1][0] + z * gray[0][1][1];
  float gray_10z = (1.F - z) * gray[1][0][0] + z * gray[1][0][1];
  float gray_11z = (1.F - z) * gray[1][1][0] + z * gray[1][1][1];
  
  float gray_0y = (1.F - y) * gray_00z + y * gray_01z;
  float gray_1y = (1.F - y) * gray_10z + y * gray_11z;
  
  interpolated_gray = (1.F - x) * gray_0y + x * gray_1y;
  
  return (unsigned short int)max(0, min(1<<16, (int)interpolated_gray));
}
*/

unsigned short int myglInterpolatedGray (int site_i, int site_j, int site_k)
{
  float x, y, z;
  float dx, dy, dz;
  float weight, sum1, sum2;
  float interpolated_gray;
  
  int voxel_i, voxel_j, voxel_k;
  int i, j, k;
  
  
  x = site_i + 0.5F;
  y = site_j + 0.5F;
  z = site_k + 0.5F;
  
  voxel_i = (int)((float)site_i * mygl.scale_inv_x);
  voxel_j = (int)((float)site_j * mygl.scale_inv_y);
  voxel_k = (int)((float)site_k * mygl.scale_inv_z);
  
  sum1 = sum2 = 0.F;
  
  for (i = max(0, voxel_i-voxel_di); i <= min(voxel_i+voxel_di, mygl.input_image_pix_x-1); i++)
    {
      dx = (float)(i + 0.5F) * mygl.scale_x - x;
      dx *= dx;
      
      for (j = max(0, voxel_j-voxel_dj); j <= min(voxel_j+voxel_dj, mygl.input_image_pix_y-1); j++)
	{
	  dy = (float)(j + 0.5F) * mygl.scale_y - y;
	  dy *= dy;
	  
	  for (k = max(0, voxel_k-voxel_dk); k <= min(voxel_k+voxel_dk, mygl.input_slices-1); k++)
	    {
	      dz = (float)(k + 0.5F) * mygl.scale_z - z;
	      dz *= dz;
	      
	      weight = dx + dy + dz;
	      
	      if (weight > res_factor_par1) continue;
	      
	      weight = 1.F - weight * res_factor_par2;
	      sum1 += weight * (float)mygl.medical_data[ myglVoxelId(i, j, k) ];
	      sum2 += weight;
	    }
	}
    }
  interpolated_gray = sum1 / sum2;
  
  return (unsigned short int)max(0, min(1<<16, (int)interpolated_gray));
}


unsigned int *myglSiteDataPointer (int site_i, int site_j, int site_k, MyGL *mygl)
{
  int i, j, k, ii, jj, kk;

  DataBlock *data_block_p;

  
  if (site_i <= -1 || site_i >= mygl->sites_x ||
      site_j <= -1 || site_j >= mygl->sites_y ||
      site_k <= -1 || site_k >= mygl->sites_z)
    {
      return NULL;
    }
  i = site_i >> mygl->shift;
  j = site_j >> mygl->shift;
  k = site_k >> mygl->shift;
  
  data_block_p = &mygl->data_block[(i * mygl->blocks_y + j) * mygl->blocks_z + k];
  
  if (data_block_p->site_data == NULL)
    {
      return NULL;   // an empty block (solid nodes region) is addressed
    }
  else
    {
      ii = site_i - (i << mygl->shift);
      jj = site_j - (j << mygl->shift);
      kk = site_k - (k << mygl->shift);
      
      return &data_block_p->site_data[(((ii << mygl->shift) + jj) << mygl->shift) + kk];
    }
}


void myglRotate (float x0, float y0, float z0,
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


void myglRotate (float x0, float y0, float z0,
		 float sn1, float cs1, float sn2, float cs2,
		 float *x1, float *y1, float *z1)
{
  float temp;
  
  
  temp = cs2 * z0 - sn2 * y0;

  *x1 = sn1 * temp + cs1 * x0;
  *y1 = sn2 * z0   + cs2 * y0;
  *z1 = cs1 * temp - sn1 * x0;
}


void myglAntiRotate (float x0, float y0, float z0,
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


void EditLastTriangle (float t_x, float t_y, float longitude, float latitude, float size_factor)
{
  float cx, cy, cz;
  float nx0, ny0, nz0;
  float nx1, ny1, nz1;
  float x, y, z;
  
  int i;
  
  Triangle *triangle_p;
  
  
  triangle_p = &mygl.boundary[ last_triangle.boundary_id ].triangle[ last_triangle.triangle_id ];
  
  cx = (1.F / 3.F) * (triangle_p->v[0].pos_x + triangle_p->v[1].pos_x + triangle_p->v[2].pos_x);
  cy = (1.F / 3.F) * (triangle_p->v[0].pos_y + triangle_p->v[1].pos_y + triangle_p->v[2].pos_y);
  cz = (1.F / 3.F) * (triangle_p->v[0].pos_z + triangle_p->v[1].pos_z + triangle_p->v[2].pos_z);
  
  for (i = 0; i < 3; i++)
    {
      nx0 = triangle_p->v[i].pos_x - cx;
      ny0 = triangle_p->v[i].pos_y - cy;      
      nz0 = triangle_p->v[i].pos_z - cz;
      
      myglRotate (nx0, ny0, nz0, longitude, latitude, &nx1, &ny1, &nz1);
      
      triangle_p->v[i].pos_x = size_factor * nx1 + cx;
      triangle_p->v[i].pos_y = size_factor * ny1 + cy;
      triangle_p->v[i].pos_z = size_factor * nz1 + cz;
      
      myglTransformVertex (triangle_p->v[i].pos_x, triangle_p->v[i].pos_y, triangle_p->v[i].pos_z, &x, &y, &z);
      
      x += t_x;
      y += t_y;
      
      myglAntiTransformVertex (x, y, z, &triangle_p->v[i].pos_x, &triangle_p->v[i].pos_y, &triangle_p->v[i].pos_z);
    }
}


void ChangeTrianglePars (int boundary_id, int triangle_id, float dp_avg, float dp_amp, float dp_phs)
{
  Triangle *triangle_p = &mygl.boundary[ boundary_id ].triangle[ triangle_id ];
  
  triangle_p->pressure_avg += dp_avg;
  triangle_p->pressure_avg = fmaxf (0.0, triangle_p->pressure_avg);
  
  triangle_p->pressure_amp += dp_amp;
  triangle_p->pressure_amp = fminf (triangle_p->pressure_avg, triangle_p->pressure_amp);
  
  triangle_p->pressure_phs += dp_phs;
  triangle_p->pressure_phs = fmaxf (0.0, triangle_p->pressure_phs);
  triangle_p->pressure_phs = fminf (180.0, triangle_p->pressure_phs);
}


void DisplayTrianglePars (int boundary_id, int triangle_id)
{
  float x1, y1, z1;
  float x2, y2;
  
  int n;
  
  char pars_string[256];
  char pressure_avg[256];
  char pressure_amp[256];
  char pressure_phs[256];
  
  Triangle *triangle_p;

  
  x2 = -1.e+30F;
  y2 = -1.e+30F;
  
  triangle_p = &mygl.boundary[ boundary_id ].triangle[ triangle_id ];
  
  for (n = 0; n < 3; n++)
    {
      myglTransformVertex (triangle_p->v[n].pos_x, triangle_p->v[n].pos_y, triangle_p->v[n].pos_z, &x1, &y1, &z1);
      
      x2 = fmaxf(x2, x1);
      y2 = fmaxf(y2, y1);
    }
  sprintf (pressure_avg, "%.0f,", triangle_p->pressure_avg);
  sprintf (pressure_amp, "%.1f,", triangle_p->pressure_amp);
  sprintf (pressure_phs, "%.0f", triangle_p->pressure_phs);
  
  strcpy (pars_string, pressure_avg);
  strcat (pars_string, pressure_amp);
  strcat (pars_string, pressure_phs);
  
  myglDisplayString ((int)x2, (int)y2, pars_string, GLUT_BITMAP_TIMES_ROMAN_24);
}


void DisplaySite (int site_i, int site_j, int site_k)
{
  float x1, y1, z1;
  float x2, y2, z2;
  
  char indices_string[256];
  char coord_i[256];
  char coord_j[256];
  char coord_k[256];
  
  
  x1 = (float)site_i * mygl.lattice_to_system - mygl.half_dim_x;
  y1 = (float)site_j * mygl.lattice_to_system - mygl.half_dim_y;
  z1 = (float)site_k * mygl.lattice_to_system - mygl.half_dim_z;
  
  myglTransformVertex (x1, y1, z1, &x2, &y2, &z2);
  
  sprintf (coord_i, "%i,", site_i);
  sprintf (coord_j, "%i,", site_j);
  sprintf (coord_k, "%i", site_k);
  
  strcat (indices_string, coord_i);
  strcat (indices_string, coord_j);
  strcat (indices_string, coord_k);
  
  myglDisplayString ((int)x2, (int)y2, indices_string, GLUT_BITMAP_TIMES_ROMAN_24);
}
  

void DeleteLastTriangle (void)
{
  int n;
  
  Boundary *boundary_p;
  
  Triangle *triangle1, *triangle2;
  
  
  boundary_p = &mygl.boundary[ last_triangle.boundary_id ];
  
  if (last_triangle.triangle_id == --boundary_p->triangles)
    {
      return;
    }
  triangle2 = &boundary_p->triangle[ last_triangle.triangle_id ];
  triangle1 = &boundary_p->triangle[ boundary_p->triangles ];
  
  for (n = 0; n < 3; n++)
    {
      triangle2->v[n].pos_x = triangle1->v[n].pos_x;
      triangle2->v[n].pos_y = triangle1->v[n].pos_y;
      triangle2->v[n].pos_z = triangle1->v[n].pos_z;
    }
  triangle2->boundary_id     = triangle1->boundary_id;
  triangle2->boundary_dir    = triangle1->boundary_dir;
  triangle2->boundary_config = triangle1->boundary_config;
  
  triangle2->pressure_avg = triangle1->pressure_avg;
  triangle2->pressure_amp = triangle1->pressure_amp;
  triangle2->pressure_phs = triangle1->pressure_phs;
}


void myglInitRayTracing (void)
{
  float block_size_inv;
  
  int triangle_index;
  int n;
  int ii1, jj1, kk1;
  int ii2, jj2, kk2;
  int ii3, jj3, kk3;
  int iia, jja, kka;
  int iib, jjb, kkb;
  int ii0, jj0, kk0;
  
  unsigned int boundary_index;
  
  Triangle *triangle_p;
  
  LinkedTriangle *new_linked_triangle;
  
  DataBlock *data_block_p;
  
  
  for (n = 0; n < mygl.blocks; n++)
    {
      mygl.data_block[ n ].linked_triangle = NULL;
    }
  
  block_size_inv = mygl.blocks_x / mygl.dim_x;
  
  for (boundary_index = 0; boundary_index < BOUNDARIES; boundary_index++)
    {
      for (triangle_index = 0; triangle_index < mygl.boundary[ boundary_index ].triangles; triangle_index++)
	{
	  triangle_p = &mygl.boundary[ boundary_index ].triangle[ triangle_index ];
	  
	  ii1 = (int)(block_size_inv * (triangle_p->v[0].pos_x + mygl.half_dim_x));
	  jj1 = (int)(block_size_inv * (triangle_p->v[0].pos_y + mygl.half_dim_y));
	  kk1 = (int)(block_size_inv * (triangle_p->v[0].pos_z + mygl.half_dim_z));
	  
	  ii2 = (int)(block_size_inv * (triangle_p->v[1].pos_x + mygl.half_dim_x));
	  jj2 = (int)(block_size_inv * (triangle_p->v[1].pos_y + mygl.half_dim_y));
	  kk2 = (int)(block_size_inv * (triangle_p->v[1].pos_z + mygl.half_dim_z));
	  
	  ii3 = (int)(block_size_inv * (triangle_p->v[2].pos_x + mygl.half_dim_x));
	  jj3 = (int)(block_size_inv * (triangle_p->v[2].pos_y + mygl.half_dim_y));
	  kk3 = (int)(block_size_inv * (triangle_p->v[2].pos_z + mygl.half_dim_z));
	  
	  iia = min(ii1, min(ii2, ii3));
	  iib = max(ii1, max(ii2, ii3));
	  jja = min(jj1, min(jj2, jj3));
	  jjb = max(jj1, max(jj2, jj3));
	  kka = min(kk1, min(kk2, kk3));
	  kkb = max(kk1, max(kk2, kk3));
	  
	  for (ii0 = iia; ii0 <= iib; ii0++)
	    {
	      if (ii0 < 0 || ii0 >= mygl.blocks_x) continue;
	      
	      for (jj0 = jja; jj0 <= jjb; jj0++)
		{
		  if (jj0 < 0 || jj0 >= mygl.blocks_y) continue;
		  
		  for (kk0 = kka; kk0 <= kkb; kk0++)
		    {
		      if (kk0 < 0 || kk0 >= mygl.blocks_z) continue;
		      
		      data_block_p = &mygl.data_block[ (ii0 * mygl.blocks_y + jj0) * mygl.blocks_z + kk0 ];
		      
		      if (data_block_p->site_data == NULL)
			{
			  continue;
			}
		      if (data_block_p->linked_triangle == NULL)
			{
			  data_block_p->linked_triangle = (LinkedTriangle *)malloc(sizeof(LinkedTriangle));
			  data_block_p->linked_triangle->boundary_id = boundary_index;
			  data_block_p->linked_triangle->triangle_id = triangle_index;
			  data_block_p->linked_triangle->next        = NULL;
			}
		      else
			{
			  new_linked_triangle = (LinkedTriangle *)malloc(sizeof(LinkedTriangle));
			  new_linked_triangle->boundary_id = boundary_index;
			  new_linked_triangle->triangle_id = triangle_index;
			  new_linked_triangle->next        = data_block_p->linked_triangle->next;
			  
			  data_block_p->linked_triangle->next = new_linked_triangle;
			}
		    }
		}
	    }
	}
    }
}


int myglIsSegmentIntercepted (int block_id1, int block_id2,
			      float x, float y, float z, float nx, float ny, float nz,
			      unsigned int *site_data_p)
{
  float t, v, w;
  float t_max;
  
  int is_intersected;
  int n;
  int block_id[2], blocks;
  
  Triangle *triangle_p;
  
  LinkedTriangle *lt_p;
  
  DataBlock *data_block_p;
  
  
  blocks = 0;
  
  if (mygl.data_block[ block_id1 ].site_data != NULL)
    {
      block_id[ 0 ] = block_id1;
      ++blocks;
    }
  if (mygl.data_block[ block_id2 ].site_data != NULL)
    {
      block_id[ blocks ] = block_id2;
      ++blocks;
    }
  
  is_intersected = 0;
  
  if (blocks == 0)
    {
      return is_intersected;
    }
  
  t_max = 1.F * mygl.lattice_to_system;
  
  for (n = 0; n < blocks; n++)
    {
      data_block_p = &mygl.data_block[ block_id[n] ];
      
      for (lt_p = data_block_p->linked_triangle;
	   lt_p != NULL;
	   lt_p = lt_p->next)
	{
	  triangle_p = &mygl.boundary[ lt_p->boundary_id ].triangle[ lt_p->triangle_id ];
	  
	  if (!myglRayTriangleIntersection (x, y, z, nx, ny, nz,
					    triangle_p->v[0].pos_x, triangle_p->v[0].pos_y, triangle_p->v[0].pos_z,
					    triangle_p->v[1].pos_x, triangle_p->v[1].pos_y, triangle_p->v[1].pos_z,
					    triangle_p->v[2].pos_x, triangle_p->v[2].pos_y, triangle_p->v[2].pos_z,
					    &t, &v, &w))
	    {
	      continue;
	    }
	  
	  if (t > t_max) continue;
	  
	  is_intersected = 1;
	  t_max = t;
	  
	  if (lt_p->boundary_id == (short int)INLET_BOUNDARY)
	    {
	      *site_data_p = INLET_TYPE;
	      *site_data_p |= (triangle_p->boundary_dir << BOUNDARY_DIR_SHIFT);
	      *site_data_p |= (triangle_p->boundary_id  << BOUNDARY_ID_SHIFT);
	    }
	  else if (lt_p->boundary_id == (short int)OUTLET_BOUNDARY)
	    {
	      *site_data_p = OUTLET_TYPE;
	      *site_data_p |= (triangle_p->boundary_dir << BOUNDARY_DIR_SHIFT);
	      *site_data_p |= (triangle_p->boundary_id  << BOUNDARY_ID_SHIFT);
	    }
	  else if (lt_p->boundary_id == (short int)WALL_BOUNDARY)
	    {
	      *site_data_p = NULL_TYPE;
	    }
	}
    }
  return is_intersected;
}
   
void myglEndRayTracing (void)
{
  int n;
  
  LinkedTriangle *old_linked_triangle, *linked_triangle_p;
  
  DataBlock *data_block_p;
  
  
  for (n = 0; n < mygl.blocks; n++)
    {
      data_block_p = &mygl.data_block[ n ];
      
      for (linked_triangle_p = data_block_p->linked_triangle;
	   linked_triangle_p != NULL;)
	{
	  old_linked_triangle = linked_triangle_p;
	  linked_triangle_p = linked_triangle_p->next;
	  free(old_linked_triangle);
	}
    }
}


void myglNearestLatticeDirection (float nx0, float ny0, float nz0, unsigned int *boundary_dir)
{
  float nx1, ny1, nz1;
  float s0, s1;
  float temp;
  
  unsigned int l;
  
  
  s0 = -1.e+30F;
  
  for (l = 0; l < 14; l++)
    {
      nx1 = (float)e_x[l];
      ny1 = (float)e_y[l];
      nz1 = (float)e_z[l];      
      temp = 1.F / sqrtf(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
      nx1 *= temp;
      ny1 *= temp;
      nz1 *= temp;
      
      s1 = fabsf(nx0 * nx1 + ny0 * ny1 + nz0 * nz1);
      
      if (s1 > s0)
	{
	  s0 = s1;
	  
	  *boundary_dir = l;
	}
    }
}


void myglSetBoundaryDirections (void)
{
  float nx, ny, nz;
  
  int l;
  
  unsigned int boundary_index;
  
  Triangle *triangle_p;
  
  
  for (boundary_index = 0; boundary_index < BOUNDARIES; boundary_index++)
    {
      for (l = 0; l < mygl.boundary[ boundary_index ].triangles; l++)
	{
	  triangle_p = &mygl.boundary[ boundary_index ].triangle[ l ];
	  
	  myglTriangleNormal (triangle_p->v[0].pos_x, triangle_p->v[0].pos_y, triangle_p->v[0].pos_z,
			      triangle_p->v[1].pos_x, triangle_p->v[1].pos_y, triangle_p->v[1].pos_z,
			      triangle_p->v[2].pos_x, triangle_p->v[2].pos_y, triangle_p->v[2].pos_z,
			      &nx, &ny, &nz);
	  
	  myglNearestLatticeDirection (nx, ny, nz, &triangle_p->boundary_dir);
	  
	  triangle_p->boundary_id = l;
	}
    }
}


int myglIsSuperficialSite (int site_i, int site_j, int site_k)
{
  int neigh_i, neigh_j, neigh_k;
  int l;
  
  unsigned int *site_data_p;
  
  
  for (l = 0; l < 14; l++)
    {
      neigh_i = site_i + e_x[ l ];
      neigh_j = site_j + e_y[ l ];
      neigh_k = site_k + e_z[ l ];
      
      site_data_p = myglSiteDataPointer (neigh_i, neigh_j, neigh_k, &mygl);
      
      if (site_data_p == NULL ||
	  *site_data_p == SOLID_TYPE ||
	  *site_data_p == NULL_TYPE)
	{
	  return 1;
	}
    }
  return 0;
}


void myglAddSuperficialSite (int iters, int site_i, int site_j, int site_k)
{
  if (mygl.superficial_sites == mygl.superficial_sites_max)
    {
      mygl.superficial_sites_max <<= 1;
      mygl.superficial_site = (SuperficialSite *)realloc(mygl.superficial_site,
							 sizeof(SuperficialSite) * mygl.superficial_sites_max);
    }
  mygl.superficial_site[ mygl.superficial_sites ].i = site_i;
  mygl.superficial_site[ mygl.superficial_sites ].j = site_j;
  mygl.superficial_site[ mygl.superficial_sites ].k = site_k;
  mygl.superficial_site[ mygl.superficial_sites ].iters = iters;
  ++mygl.superficial_sites;
}


int myglInitialiseBlock (int site_i, int site_j, int site_k, int selected_gray, int iters)
{
  int block_id, site_id;
  int i, j, k;
  int ii, jj, kk;
  int m;
  
  unsigned int *site_data_p;
  
  DataBlock *data_block_p;
  
  
  i = site_i >> mygl.shift;
  j = site_j >> mygl.shift;
  k = site_k >> mygl.shift;
  block_id = (i * mygl.blocks_y + j) * mygl.blocks_z + k;
  
  i <<= mygl.shift;
  j <<= mygl.shift;
  k <<= mygl.shift;
  
  data_block_p = &mygl.data_block[ block_id ];
  
  if (data_block_p->site_data == NULL)
    {
      data_block_p->site_data = (unsigned int *)malloc(sizeof(unsigned int) * mygl.sites_in_a_block);
      
      data_block_p->site_iters = (unsigned short int *)malloc(sizeof(unsigned short int) * mygl.sites_in_a_block);
      
      m = -1;
      
      for (ii = i; ii < i + mygl.block_size; ii++)
	{
	  for (jj = j; jj < j + mygl.block_size; jj++)
	    {
	      for (kk = k; kk < k + mygl.block_size; kk++)
		{
		  if (myglInterpolatedGray (ii, jj, kk) > selected_gray)
		    {
		      data_block_p->site_data[ ++m ] = NULL_TYPE;
		    }
		  else
		    {
		      data_block_p->site_data[ ++m ] = SOLID_TYPE;
		    }
		}
	    }
	}
      data_block_p->is_void = 1;
    }
  ii = site_i - i;
  jj = site_j - j;
  kk = site_k - k;
  site_id = (((ii << mygl.shift) + jj) << mygl.shift) + kk;
  
  site_data_p = &data_block_p->site_data[ site_id ];
  
  if (*site_data_p != NULL_TYPE)
    {
      return 0;
    }
  *site_data_p = FLUID_TYPE;
  data_block_p->site_iters[ site_id ] = iters;
  
  data_block_p->is_void = 0;
  
  return 1;
}


void myglDeallocateBlock (int n)
{
  DataBlock *data_block_p = &mygl.data_block[ n ];
  
  if (data_block_p->site_data != NULL)
    {
      free(data_block_p->site_data);
      data_block_p->site_data = NULL;
    }
  if (data_block_p->site_iters != NULL)
    {
      free(data_block_p->site_iters);
      data_block_p->site_iters = NULL;
    }
  data_block_p->site_label = NULL;
  data_block_p->is_void = 1;
}


void myglReconstructSystem (int selected_pixel_x,
			    int selected_pixel_y,
			    int selected_slice,
			    unsigned short int selected_gray)
{
  float seconds;
  
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int l, n;
  int iterations;
  int sites_a, sites_b;
  int index_a;
  int are_fluid_sites_incrementing;
  int fluid_sites;
  
  SiteLocation *site_location_a_p, *site_location_b_p;
  
  
  seconds = myClock ();
  
  for (n = 0; n < mygl.blocks; n++)
    {
      myglDeallocateBlock (n);
    }
  mygl.superficial_sites = 0;
  iterations = 0;
  
  myglFromVoxelToSiteCoords (selected_pixel_x, selected_pixel_y, selected_slice, &i, &j, &k);
  
  res_factor_par1 = 4.F * res_factor * res_factor;
  res_factor_par2 = 1.F / res_factor_par1;
  
  if (!myglInitialiseBlock (i, j, k, selected_gray, iterations))
    {
      return;
    }
  site_location_a_p = &site_location_a[ 0 ];
  site_location_a_p->i = i;
  site_location_a_p->j = j;
  site_location_a_p->k = k;
  sites_a = 1;
  
  fluid_sites = 1;
  are_fluid_sites_incrementing = 1;
  
  while (are_fluid_sites_incrementing)
    {
      ++iterations;
      sites_b = 0;
      are_fluid_sites_incrementing = 0;
      
      for (index_a = 0; index_a < sites_a; index_a++)
	{
	  site_location_a_p = &site_location_a[ index_a ];
	  i = site_location_a_p->i;
	  j = site_location_a_p->j;
	  k = site_location_a_p->k;
	  
	  for (l = 0; l < 14; l++)
	    {
	      neigh_i = i + e_x[ l ];
	      neigh_j = j + e_y[ l ];
	      neigh_k = k + e_z[ l ];
	      
	      if (neigh_i == -1 || neigh_i >= mygl.sites_x) continue;
	      if (neigh_j == -1 || neigh_j >= mygl.sites_y) continue;
	      if (neigh_k == -1 || neigh_k >= mygl.sites_z) continue;
	      
	      if (!myglInitialiseBlock (neigh_i, neigh_j, neigh_k, selected_gray, iterations))
		{
		  continue;
		}
	      are_fluid_sites_incrementing = 1;
	      
	      if (sites_b == BUFFERS_SIZE)
		{
		  printf (" buffer size too small\n");
		  printf (" the execution is terminated\n");
		  exit(1);
		}
	      ++fluid_sites;
	      
	      site_location_b_p = &site_location_b[ sites_b ];
	      site_location_b_p->i = neigh_i;
	      site_location_b_p->j = neigh_j;
	      site_location_b_p->k = neigh_k;
	      ++sites_b;
	    }
	}
      --iterations;
      
      for (index_a = 0; index_a < sites_a; index_a++)
	{
	  if (myglIsSuperficialSite (site_location_a[ index_a ].i,
				     site_location_a[ index_a ].j,
				     site_location_a[ index_a ].k))
	    {
	      myglAddSuperficialSite (iterations,
				      site_location_a[ index_a ].i,
				      site_location_a[ index_a ].j,
				      site_location_a[ index_a ].k);
	    }
	}
      ++iterations;
      
      if (are_fluid_sites_incrementing)
	{
	  site_location_a_p = site_location_a;
	  site_location_a = site_location_b;
	  site_location_b = site_location_a_p;
	  sites_a = sites_b;
	}
    }
  --iterations;
  
  for (index_a = 0; index_a < sites_a; index_a++)
    {
      if (myglIsSuperficialSite (site_location_a[ index_a ].i,
				 site_location_a[ index_a ].j,
				 site_location_a[ index_a ].k))
	{
	  myglAddSuperficialSite (iterations,
				  site_location_a[ index_a ].i,
				  site_location_a[ index_a ].j,
				  site_location_a[ index_a ].k);
	}
    }
  printf ("Reconstruction, gray: %i, sites: %i, superficial sites: %i, time: %.3f\n",
	  selected_gray, fluid_sites, mygl.superficial_sites, myClock () - seconds);
  fflush (stdout);
}


void myglFluidSitesIterativeSearching (int selected_pixel_x,
				       int selected_pixel_y,
				       int selected_slice,
				       unsigned short int selected_gray)
{
  float seconds;
  float x, y, z;
  float nx, ny, nz;
  
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int l, m, n;
  int bi, bj, bk;
  int si, sj, sk;
  int block_id1, block_id2;
  int site_id;
  int iterations;
  int sites_a, sites_b;
  int index_a;
  int are_fluid_sites_incrementing;
  int fluid_sites;
  
  unsigned int *site_data_p;

  DataBlock *data_block_p;
  
  SiteLocation *site_location_a_p, *site_location_b_p;
  
  
  seconds = myClock ();
  
  for (n = 0; n < mygl.blocks; n++)
    {
      mygl.data_block[ n ].is_void = 1;
      
      if (mygl.data_block[ n ].site_data == NULL) continue;
      
      for (m = 0; m < mygl.sites_in_a_block; m++)
	{
	  if (mygl.data_block[ n ].site_data[ m ] != SOLID_TYPE)
	    {
	      mygl.data_block[ n ].site_data[ m ] = NULL_TYPE;
	    }
	}
    }
  myglInitRayTracing ();
  
  myglSetBoundaryDirections ();
  
  mygl.superficial_sites = 0;
  
  myglFromVoxelToSiteCoords (selected_pixel_x, selected_pixel_y, selected_slice, &i, &j, &k);
  
  myglFromSiteToGridCoords (i, j, k, &block_id1, &site_id);
  
  data_block_p = &mygl.data_block[ block_id1 ];
  
  data_block_p->site_data[ site_id ] = FLUID_TYPE;
  data_block_p->is_void = 0;
  
  site_location_a_p = &site_location_a[ 0 ];
  site_location_a_p->i = i;
  site_location_a_p->j = j;
  site_location_a_p->k = k;
  sites_a = 1;
  
  iterations = 0;
  data_block_p->site_iters[ site_id ] = iterations;
  
  fluid_sites = 1;
  are_fluid_sites_incrementing = 1;
  
  while (are_fluid_sites_incrementing)
    {
      ++iterations;
      
      sites_b = 0;
      
      are_fluid_sites_incrementing = 0;
      
      for (index_a = 0; index_a < sites_a; index_a++)
	{
	  site_location_a_p = &site_location_a[ index_a ];
	  i = site_location_a_p->i;
	  j = site_location_a_p->j;
	  k = site_location_a_p->k;
	  
	  bi = i >> mygl.shift;
	  bj = j >> mygl.shift;
	  bk = k >> mygl.shift;
	  block_id1 = (bi * mygl.blocks_y + bj) * mygl.blocks_z + bk;
	  
	  x = (i + 0.0F) * mygl.lattice_to_system - mygl.half_dim_x;
	  y = (j + 0.0F) * mygl.lattice_to_system - mygl.half_dim_y;
	  z = (k + 0.0F) * mygl.lattice_to_system - mygl.half_dim_z;
	  
	  for (l = 1; l < 14; l++)
	    {
	      neigh_i = i + e_x[ l ];
	      neigh_j = j + e_y[ l ];
	      neigh_k = k + e_z[ l ];
	      
	      if (neigh_i == -1 || neigh_i >= mygl.sites_x) continue;
	      if (neigh_j == -1 || neigh_j >= mygl.sites_y) continue;
	      if (neigh_k == -1 || neigh_k >= mygl.sites_z) continue;
	      
	      bi = neigh_i >> mygl.shift;
	      bj = neigh_j >> mygl.shift;
	      bk = neigh_k >> mygl.shift;
	      block_id2 = (bi * mygl.blocks_y + bj) * mygl.blocks_z + bk;
	      
	      if ((data_block_p = &mygl.data_block[ block_id2 ])->site_data == NULL)
		{
		  continue;
		}
	      si = neigh_i - (bi << mygl.shift);
	      sj = neigh_j - (bj << mygl.shift);
	      sk = neigh_k - (bk << mygl.shift);
	      site_id = (((si << mygl.shift) + sj) << mygl.shift) + sk;
	      
	      if (*(site_data_p = &data_block_p->site_data[ site_id ]) != NULL_TYPE)
		{
		  continue;
		}
	      nx = (float)e_x[ l ];
	      ny = (float)e_y[ l ];
	      nz = (float)e_z[ l ];
	      
	      if (myglIsSegmentIntercepted (block_id1, block_id2, x, y, z, nx, ny, nz, site_data_p))
		{
		  continue;
		}
	      *site_data_p = FLUID_TYPE;
	      data_block_p->site_iters[ site_id ] = iterations;
	      data_block_p->is_void = 0;
	      
	      are_fluid_sites_incrementing = 1;
	      
	      if (sites_b == BUFFERS_SIZE)
		{
		  printf (" buffer size too small\n");
		  printf (" the execution is terminated\n");
		  exit(1);
		}
	      ++fluid_sites;
	      
	      site_location_b_p = &site_location_b[ sites_b ];
	      site_location_b_p->i = neigh_i;
	      site_location_b_p->j = neigh_j;
	      site_location_b_p->k = neigh_k;
	      ++sites_b;
	    }
	}
      --iterations;
      
      for (index_a = 0; index_a < sites_a; index_a++)
	{
	  if (myglIsSuperficialSite (site_location_a[ index_a ].i,
				     site_location_a[ index_a ].j,
				     site_location_a[ index_a ].k))
	    {
	      myglAddSuperficialSite (iterations,
				      site_location_a[ index_a ].i,
				      site_location_a[ index_a ].j,
				      site_location_a[ index_a ].k);
	    }
	}
      ++iterations;
      
      if (are_fluid_sites_incrementing)
	{
	  site_location_a_p = site_location_a;
	  site_location_a = site_location_b;
	  site_location_b = site_location_a_p;
	  sites_a = sites_b;
	}
    }
  --iterations;
  
  for (index_a = 0; index_a < sites_a; index_a++)
    {
      if (myglIsSuperficialSite (site_location_a[ index_a ].i,
				 site_location_a[ index_a ].j,
				 site_location_a[ index_a ].k))
	{
	  myglAddSuperficialSite (iterations,
				  site_location_a[ index_a ].i,
				  site_location_a[ index_a ].j,
				  site_location_a[ index_a ].k);
	}
    }
  
  myglEndRayTracing ();
  
  printf ("Iterative searching, gray: %i, sites: %i, superficial sites: %i, time: %.3f\n",
	  mygl.selected_gray, fluid_sites, mygl.superficial_sites, myClock () - seconds);
  fflush (stdout);
}


void myglSetBoundaryConfigurations (void)
{
#define C 0.57735027F
  
  float seconds;
  float norm_ex[] = { 1.F,-1.F, 0.F, 0.F, 0.F, 0.F, C,-C, C,-C, C,-C, C,-C};
  float norm_ey[] = { 0.F, 0.F, 1.F,-1.F, 0.F, 0.F, C,-C, C,-C,-C, C,-C, C};
  float norm_ez[] = { 0.F, 0.F, 0.F, 0.F, 1.F,-1.F, C,-C,-C, C, C,-C,-C, C};
  float dir_x, dir_y, dir_z;
  float score, scalar_product;
  
  int unknown[13], unknowns;
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int l, m, n;
  int u;
  
  unsigned int *site_data_p;
  unsigned int *neigh_site_data_p;
  unsigned int boundary_config;
  unsigned int boundary_dir;
  
  DataBlock *data_block_p;
  
  
  seconds = myClock ();
  
  n = -1;
  
  for (i = 0; i < mygl.blocks_x; i++)
    {
      for (j = 0; j < mygl.blocks_y; j++)
	{
	  for (k = 0; k < mygl.blocks_z; k++)
	    {
	      if ((data_block_p = &mygl.data_block[ ++n ])->is_void) continue;
	      
	      m = -1;
	      
	      for (site_i = i * mygl.block_size; site_i < i * mygl.block_size + mygl.block_size; site_i++)
		{
		  for (site_j = j * mygl.block_size; site_j < j * mygl.block_size + mygl.block_size; site_j++)
		    {
		      for (site_k = k * mygl.block_size; site_k < k * mygl.block_size + mygl.block_size; site_k++)
			{
			  site_data_p = &data_block_p->site_data[ ++m ];
			  
			  if (*site_data_p == SOLID_TYPE || *site_data_p == NULL_TYPE)
			    {
			      continue;
			    }
			  boundary_config = 0U;
			  unknowns = 0;
			  
			  for (l = 0; l < 14; l++)
			    {
			      neigh_i = site_i + e_x[ l ];
			      neigh_j = site_j + e_y[ l ];
			      neigh_k = site_k + e_z[ l ];
			      
			      neigh_site_data_p = myglSiteDataPointer (neigh_i, neigh_j, neigh_k, &mygl);
			      
			      if (!(neigh_site_data_p == NULL ||
				    *neigh_site_data_p == SOLID_TYPE ||
				    *neigh_site_data_p == NULL_TYPE))
				{
				  boundary_config |= (1U << l);
				}
			      else
				{
				  unknown[ unknowns ] = inv_dir[ l ];
				  ++unknowns;
				}
			    }
			  if ((*site_data_p & SITE_TYPE_MASK) == FLUID_TYPE)
			    {
			      if (unknowns == 0) continue;
			      
			      *site_data_p |= (boundary_config << BOUNDARY_CONFIG_SHIFT);
			      
			      dir_x = 0.F;
			      dir_y = 0.F;
			      dir_z = 0.F;
			      
			      for (l = 0; l < unknowns; l++)
				{
				  u = unknown[ l ];
				  dir_x += norm_ex[ u ];
				  dir_y += norm_ey[ u ];
				  dir_z += norm_ez[ u ];
				}
			      boundary_dir = 0;
			      
			      score = 0.F;
			      
			      for (l = 0; l < unknowns; l++)
				{
				  u = unknown[ l ];
				  scalar_product = dir_x * norm_ex[ u ] + dir_y * norm_ey[ u ] + dir_z * norm_ez[ u ];
				  
				  if (scalar_product > score)
				    {
				      score        = scalar_product;
				      boundary_dir = u;
				    }
				}
			      *site_data_p |= (boundary_dir << BOUNDARY_DIR_SHIFT);
			    }
			  else
			    {
			      *site_data_p |= (boundary_config << BOUNDARY_CONFIG_SHIFT);
			      
			      boundary_dir = (*site_data_p & BOUNDARY_DIR_MASK) >> BOUNDARY_DIR_SHIFT;
			      
			      if ((boundary_config & (1U << inv_dir[ boundary_dir ])))
				{
				  *site_data_p &= ~BOUNDARY_DIR_MASK;
				  
				  *site_data_p |= (inv_dir[ boundary_dir ] << BOUNDARY_DIR_SHIFT);
				}
			    }
			}
		    }
		}
	    }
	}
    }
  printf ("Boundary configurations setup, time: %.3f\n", myClock () - seconds);
  fflush (stdout);
}


void myglRescaleSystemResolution ()
{
  int i;
  
  for (i = 0; i < mygl.blocks; i++)
    {
      myglDeallocateBlock (i);
    }
  
  
  mygl.output_image_pix_x = (int)(mygl.input_image_pix_x * res_factor);
  mygl.output_image_pix_y = (int)(mygl.input_image_pix_y * res_factor);
  
  mygl.output_slices = (int)(mygl.input_slices * (slice_size / pixel_size) *
			     (float)mygl.output_image_pix_x / (float)mygl.input_image_pix_x);
  
  mygl.scale_x = (float)mygl.output_image_pix_x / (float)mygl.input_image_pix_x;
  mygl.scale_y = (float)mygl.output_image_pix_y / (float)mygl.input_image_pix_y;
  mygl.scale_z = (float)mygl.output_slices / (float)mygl.input_slices;
  
  mygl.scale_inv_x = 1.F / mygl.scale_x;
  mygl.scale_inv_y = 1.F / mygl.scale_y;
  mygl.scale_inv_z = 1.F / mygl.scale_z;
  
  voxel_di = 2;
  voxel_dj = 2;
  voxel_dk = nint(2.F * pixel_size / slice_size);
  
  printf (" output_slices: %i\n", mygl.output_slices);
  fflush (stdout);
  
  // system parameters setup
  
  mygl.lattice_to_system = 1.F;
  
  mygl.block_size = 8;
  mygl.shift = 3;
  
  mygl.blocks_x = mygl.output_image_pix_x >> mygl.shift;
  mygl.blocks_y = mygl.output_image_pix_y >> mygl.shift;
  mygl.blocks_z = mygl.output_slices      >> mygl.shift;
  
  if ((mygl.blocks_x << mygl.shift) < mygl.output_image_pix_x) ++mygl.blocks_x;
  if ((mygl.blocks_y << mygl.shift) < mygl.output_image_pix_y) ++mygl.blocks_y;
  if ((mygl.blocks_z << mygl.shift) < mygl.output_slices     ) ++mygl.blocks_z;
  
  mygl.sites_x = mygl.blocks_x * mygl.block_size;
  mygl.sites_y = mygl.blocks_y * mygl.block_size;
  mygl.sites_z = mygl.blocks_z * mygl.block_size;
  
  mygl.dim_x = (float)mygl.sites_x * mygl.lattice_to_system;
  mygl.dim_y = (float)mygl.sites_y * mygl.lattice_to_system;
  mygl.dim_z = (float)mygl.sites_z * mygl.lattice_to_system;
  
  mygl.half_dim_x = 0.5F * mygl.dim_x;
  mygl.half_dim_y = 0.5F * mygl.dim_y;
  mygl.half_dim_z = 0.5F * mygl.dim_z;
  
  mygl.system_size = fmaxf(mygl.dim_x, fmaxf(mygl.dim_y, mygl.dim_z));
  
  mygl.sites_in_a_block = mygl.block_size * mygl.block_size * mygl.block_size;
  
  mygl.blocks = mygl.blocks_x * mygl.blocks_y * mygl.blocks_z;
  
  
  mygl.data_block = (DataBlock *)realloc(mygl.data_block, sizeof(DataBlock) * mygl.blocks);
  
  for (i = 0; i < mygl.blocks; i++)
    {
      mygl.data_block[ i ].site_data = NULL;
      mygl.data_block[ i ].site_label = NULL;
    }
}


void myglRescaleViewpoint (float scale)
{
  ortho_x *= scale;
  ortho_y *= scale;
  viewpoint_radius *= scale;
  
  Projection ();
}


void myglRescaleTriangles (float scale)
{
  for (int n = 0; n < BOUNDARIES; n++)
    {
      for (int m = 0; m < mygl.boundary[ n ].triangles; m++)
	{
	  for (int l = 0; l < 3; l++)
	    {
	      mygl.boundary[ n ].triangle[ m ].v[l].pos_x *= scale;
	      mygl.boundary[ n ].triangle[ m ].v[l].pos_y *= scale;
	      mygl.boundary[ n ].triangle[ m ].v[l].pos_z *= scale;
	    }
	}
    }
}


void EstimateBoundaryNormal (int site_i, int site_j, int site_k, float nor[], float *length)
{
  float segment_length = 0.25F;
  float org[3];
  float x, y, z;
  float nx, ny, nz;
  float temp;
  
  int neigh_i, neigh_j, neigh_k;
  int l, n;
  int bi, bj, bk;
  int si, sj, sk;
  int block_id, site_id;
  int sites_a;
  int index_a;
  int are_fluid_sites_incrementing;
  int iterations;
  
  short int *site_label_p;
  short int global_counter;
  
  unsigned int *site_data_p, *neigh_site_data_p;
  
  DataBlock *data_block_p;
  
  SiteLocation *site_location_a_p;
  
  
  org[0] = (float)site_i;
  org[1] = (float)site_j;
  org[2] = (float)site_k;
  
  nor[0] = 0.F;
  nor[1] = 0.F;
  nor[2] = 0.F;
  
  global_counter = 1;
  mygl.stored_blocks = 0;
  
  site_location_a_p = &site_location_a[ 0 ];
  site_location_a_p->i = site_i;
  site_location_a_p->j = site_j;
  site_location_a_p->k = site_k;
  sites_a = 1;
  
  bi = site_i >> mygl.shift;
  bj = site_j >> mygl.shift;
  bk = site_k >> mygl.shift;
  block_id = (bi * mygl.blocks_y + bj) * mygl.blocks_z + bk;
  
  data_block_p = &mygl.data_block[ block_id ];
  data_block_p->site_label = (short int *)malloc(sizeof(short int) * mygl.sites_in_a_block);
  
  for (n = 0; n < mygl.sites_in_a_block; n++)
    {
      data_block_p->site_label[ n ] = 0;
    }
  si = site_i - (bi << mygl.shift);
  sj = site_j - (bj << mygl.shift);
  sk = site_k - (bk << mygl.shift);
  site_id = (((si << mygl.shift) + sj) << mygl.shift) + sk;
  
  data_block_p->site_label[ site_id ] = global_counter;
  
  are_fluid_sites_incrementing = 1;
  iterations = 0;
  
  while (++iterations <= 10 && are_fluid_sites_incrementing)
    {
      are_fluid_sites_incrementing = 0;
      
      for (index_a = 0; index_a < sites_a; index_a++)
	{
	  site_location_a_p = &site_location_a[ index_a ];
	  
	  for (l = 0; l < 14; l++)
	    {
	      neigh_i = site_location_a_p->i + e_x[ l ];
	      neigh_j = site_location_a_p->j + e_y[ l ];
	      neigh_k = site_location_a_p->k + e_z[ l ];
	      
	      if (neigh_i == -1 || neigh_i >= mygl.sites_x) continue;
	      if (neigh_j == -1 || neigh_j >= mygl.sites_y) continue;
	      if (neigh_k == -1 || neigh_k >= mygl.sites_z) continue;
	      
	      bi = neigh_i >> mygl.shift;
	      bj = neigh_j >> mygl.shift;
	      bk = neigh_k >> mygl.shift;
	      block_id = (bi * mygl.blocks_y + bj) * mygl.blocks_z + bk;
	      
	      if ((data_block_p = &mygl.data_block[ block_id ])->site_data == NULL)
		{
		  continue;
		}
	      si = neigh_i - (bi << mygl.shift);
	      sj = neigh_j - (bj << mygl.shift);
	      sk = neigh_k - (bk << mygl.shift);
	      site_id = (((si << mygl.shift) + sj) << mygl.shift) + sk;
	      
	      site_data_p = &data_block_p->site_data[ site_id ];
	      
	      if (*site_data_p == SOLID_TYPE || *site_data_p == NULL_TYPE)
		{
		  continue;
		}
	      if (!myglIsSuperficialSite (neigh_i, neigh_j, neigh_k))
		{
		  continue;
		}
	      if (data_block_p->site_label == NULL)
		{
		  data_block_p->site_label = (short int *)malloc(sizeof(short int) * mygl.sites_in_a_block);
		  
		  for (n = 0; n < mygl.sites_in_a_block; n++)
		    {
		      data_block_p->site_label[ n ] = 0;
		    }
		  if (mygl.stored_blocks == mygl.stored_blocks_max)
		    {
		      mygl.stored_blocks_max <<= 1;
		      mygl.stored_block = (int *)realloc(mygl.stored_block, sizeof(int) * mygl.stored_blocks_max);
		    }
		  mygl.stored_block[ mygl.stored_blocks ] = block_id;
		  ++mygl.stored_blocks;
		}
	      site_label_p = &data_block_p->site_label[ site_id ];
	      
	      if (*site_label_p == global_counter) continue;
	      
	      *site_label_p = global_counter;
	      
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
    }
  for (n = 0; n < mygl.stored_blocks; n++)
    {
      block_id = mygl.stored_block[ n ];
      
      free(mygl.data_block[ block_id ].site_label);
    }
  
  temp = 1.F / sqrtf(nor[0] * nor[0] + nor[1] * nor[1] + nor[2] * nor[2]);
  nor[0] *= temp;
  nor[1] *= temp;
  nor[2] *= temp;
  
  x = org[0] + 2.F * nor[0];
  y = org[1] + 2.F * nor[1];
  z = org[2] + 2.F * nor[2];
  
  neigh_i = (int)x;
  neigh_j = (int)y;
  neigh_k = (int)z;
  neigh_site_data_p = myglSiteDataPointer (neigh_i, neigh_j, neigh_k, &mygl);
  
  if (neigh_site_data_p == NULL ||
      *neigh_site_data_p == SOLID_TYPE ||
      *neigh_site_data_p == NULL_TYPE)
    {
      nor[0] = -nor[0];
      nor[1] = -nor[1];
      nor[2] = -nor[2];
    }
  nor[0] *= segment_length;
  nor[1] *= segment_length;
  nor[2] *= segment_length;
  
  x = org[0];
  y = org[1];
  z = org[2];
  
  iterations = 0;
  
  site_data_p = myglSiteDataPointer (site_i, site_j, site_k, &mygl);
  
  while (site_data_p != NULL &&
	 *site_data_p != SOLID_TYPE &&
	 *site_data_p != NULL_TYPE)
    {
      ++iterations;
      
      x += nor[0];
      y += nor[1];
      z += nor[2];
      
      neigh_i = (int)x;
      neigh_j = (int)y;
      neigh_k = (int)z;
      
      site_data_p = myglSiteDataPointer (neigh_i, neigh_j, neigh_k, &mygl);
    }
  --iterations;
  
  if (iterations == 0)
    {
      *length = -1.F;
    }
  else
    {
      *length = (float)iterations;
    }
}


int myglCreateOptimizedTriangle (unsigned int site_type,
				 int site_i, int site_j, int site_k)
{
  float triangle_factor = 1.F;
  float org[3], nor[3];
  float x, y, z;
  float longitude, latitude;
  float triangle_size;
  float length;
  
  int neigh_i, neigh_j, neigh_k;
  int l, m, n;
  int iterations, iterations_max;
  int longitude_id, latitude_id;
  int triangle_id;
  
  unsigned int *site_data_p, *neigh_site_data_p;
  
  Triangle *triangle_p;
  
  
  if (mygl.boundary[ site_type ].triangles == (int)(1U << BOUNDARY_ID_BITS))
    {
      return -1;
    }
  triangle_id = mygl.boundary[ site_type ].triangles;
  triangle_p = &mygl.boundary[ site_type ].triangle[ triangle_id ];
  
  org[0] = (float)site_i;
  org[1] = (float)site_j;
  org[2] = (float)site_k;
  
  EstimateBoundaryNormal (site_i, site_j, site_k, nor, &length);
  
  if (length < 0.F)
    {
      return -1;
    }
  
  triangle_size = triangle_factor * length * mygl.lattice_to_system;
  
  // segment center
  org[0] = 0.5F * (org[0] + (org[0] + length * nor[0]));
  org[1] = 0.5F * (org[1] + (org[1] + length * nor[1]));
  org[2] = 0.5F * (org[2] + (org[2] + length * nor[2]));
  
  site_i = (int)org[0];
  site_j = (int)org[1];
  site_k = (int)org[2];
  
  site_data_p = myglSiteDataPointer (site_i, site_j, site_k, &mygl);
  
  iterations_max = 0;
  longitude_id = 0;
  latitude_id = 0;
  
  latitude = 0.F;
  
  for (n = 0; n < 180; n++)
    {
      longitude = 0.F;
      
      for (m = 0; m < 360; m++)
	{
	  myglRotate (0.F, 0.F, 1.F,
		      longitude, latitude,
		      &nor[0], &nor[1], &nor[2]);
	  
	  x = org[0];
	  y = org[1];
	  z = org[2];
	  
	  iterations = 0;
	  
	  neigh_site_data_p = site_data_p;
	  
	  while (neigh_site_data_p != NULL &&
		 *neigh_site_data_p != SOLID_TYPE &&
		 *neigh_site_data_p != NULL_TYPE)
	    {
	      ++iterations;
	  
	      x += nor[0];
	      y += nor[1];
	      z += nor[2];
	      
	      neigh_i = (int)x;
	      neigh_j = (int)y;
	      neigh_k = (int)z;
	      
	      neigh_site_data_p = myglSiteDataPointer (neigh_i, neigh_j, neigh_k, &mygl);
	    }
	  --iterations;
	  
	  if (iterations > iterations_max)
	    {
	      iterations_max = iterations;
	      longitude_id = m;
	      latitude_id = n;
	    }
	  longitude += DEG_TO_RAD;
	}
      latitude += 2.F * DEG_TO_RAD;
    }
  if (iterations_max == 0)
    {
      return -1;
    }
  
  longitude = longitude_id * DEG_TO_RAD;
  latitude = latitude_id * 2.F * DEG_TO_RAD;
  
  myglRotate (0.F, 0.F, 1.F,
	      longitude, latitude,
	      &nor[0], &nor[1], &nor[2]);
  
  org[0] = org[0] * mygl.lattice_to_system - mygl.half_dim_x;
  org[1] = org[1] * mygl.lattice_to_system - mygl.half_dim_y;
  org[2] = org[2] * mygl.lattice_to_system - mygl.half_dim_z;
  
  myglRotate (0.F, triangle_size, 0.F,
	      longitude, latitude,
	      &triangle_p->v[0].pos_x, &triangle_p->v[0].pos_y, &triangle_p->v[0].pos_z);
  
  myglRotate (-(triangle_size / sqrtf(2.F)), -(triangle_size / sqrtf(2.F)), 0.F,
	      longitude, latitude,
	      &triangle_p->v[1].pos_x, &triangle_p->v[1].pos_y, &triangle_p->v[1].pos_z);
  
  myglRotate (+(triangle_size / sqrtf(2.F)), -(triangle_size / sqrtf(2.F)), 0.F,
	      longitude, latitude,
	      &triangle_p->v[2].pos_x, &triangle_p->v[2].pos_y, &triangle_p->v[2].pos_z);
  
  for (l = 0; l < 3; l++)
    {
      triangle_p->v[l].pos_x += org[0];
      triangle_p->v[l].pos_y += org[1];
      triangle_p->v[l].pos_z += org[2];
    }
  triangle_p->pressure_avg = 80.0F;
  triangle_p->pressure_amp = 0.0F;
  triangle_p->pressure_phs = 0.0F;
  
  ++mygl.boundary[ site_type ].triangles;
  
  return triangle_id;
}


void DisplayNull (void)
{
  ;
}

void DisplaySlice (void)
{
  float x, y;
  float delta_x, delta_y;
  float gray;
  
  int i, j, k;
  
  
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  delta_x = 2.F * screen.max_x / (float)mygl.input_image_pix_x;
  delta_y = 2.F * screen.max_y / (float)mygl.input_image_pix_y;
  
  k = mygl.selected_slice;
  
  x = -screen.max_x + 0.5F * delta_x;
  
  for (i = 0; i < mygl.input_image_pix_x - 1; i++)
    {
      glBegin (GL_TRIANGLE_STRIP);
      
      y = -screen.max_y + 0.5F * delta_y;
      
      for (j = 0; j < mygl.input_image_pix_y; j++)
	{
	  gray = mygl.medical_data[ myglVoxelId(i,j,k) ] * (1.F / (1 << 8));
	  
	  glColor3f (gray, gray, gray);	  
  	  glVertex3f (x, y, 0.F);
  	  
	  gray = mygl.medical_data[ myglVoxelId(i+1,j,k) ] * (1.F / (1 << 8));
	  
  	  glColor3f (gray, gray, gray);
  	  glVertex3f (x + delta_x, y, 0.F);
	  
  	  y += delta_y;
  	}
      x += delta_x;
      
      glEnd ();
    }
  
  glutSwapBuffers ();
}


void DisplaySystemSlow (void)
{
  // slow version!
  
  float block_size;
  float x0, y0, z0;
  float x1, y1, z1;
  float x2, y2, z2;
  
  int pixel_i, pixel_j;
  int i, j, k;
  int ii, jj, kk;
  int m, n;
  int triangle_index;
  
  unsigned int boundary_index;
  
  unsigned int site_data;

  DataBlock *data_block_p;
  
  Triangle *triangle_p;
  
  ScreenVoxel *screen_voxel_p;
  
  
  //display_id = 0;
  
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  
  glColor3f (0.F, 0.F, 0.F);
  
  glLineWidth (3.F);
  glBegin (GL_LINE_LOOP);
  
  myglTransformVertex (-mygl.half_dim_x, -mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (-mygl.half_dim_x, -mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, -mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, -mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  
  glEnd ();
  glBegin (GL_LINE_LOOP);

  myglTransformVertex (-mygl.half_dim_x, +mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (-mygl.half_dim_x, +mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, +mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, +mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  
  glEnd ();
  
  
  block_size = (float)mygl.block_size * mygl.lattice_to_system;
  
  glPointSize ((float)point_size);
  glBegin (GL_POINTS);
  
  n = -1;
  
  x0 = -mygl.half_dim_x;

  for (i = 0; i < mygl.blocks_x; i++)
    {
      y0 = -mygl.half_dim_y;
      
      for (j = 0; j < mygl.blocks_y; j++)
	{
	  z0 = -mygl.half_dim_z;
	  
	  for (k = 0; k < mygl.blocks_z; k++)
	    {
	      data_block_p = &mygl.data_block[ ++n ];
	      
	      if (data_block_p->site_data == NULL ||
		  data_block_p->is_void)
		{
		  z0 += block_size;
		  continue;
		}
	      m = -1;
	      
	      x1 = x0;
	      
	      for (ii = 0; ii < mygl.block_size; ii++)
		{
		  if ((i << mygl.shift) + ii >= mygl.sites_x) break;
	  
		  y1 = y0;
		  
		  for (jj = 0; jj < mygl.block_size; jj++)
		    {
		      if ((j << mygl.shift) + jj >= mygl.sites_y) break;
		      
		      z1 = z0;
		      
		      for (kk = 0; kk < mygl.block_size; kk++)
			{
			  if ((k << mygl.shift) + kk >= mygl.sites_z) break;
			  
			  site_data = data_block_p->site_data[ ++m ];
			  
			  if (site_data == SOLID_TYPE ||
			      site_data == NULL_TYPE)
			    {
			      z1 += mygl.lattice_to_system;
			      continue;
			    }
			  site_data = (site_data & SITE_TYPE_MASK);
			  
			  if (site_data == FLUID_TYPE)
			    {
			      glColor3f (RedComponent (data_block_p->site_iters[m]), 0.F, 0.F);
			    }
			  else if (site_data == INLET_TYPE)
			    {
			      glColor3f (0.F, 0.F, 0.F);
			    }
			  else if (site_data == OUTLET_TYPE)
			    {
			      glColor3f (0.F, 1.F, 0.F);
			    }
			  myglTransformVertex (x1, y1, z1, &x2, &y2, &z2);
			  
			  glVertex3f (x2, y2, z2);
			  
			  z1 += mygl.lattice_to_system;
			}
		      y1 += mygl.lattice_to_system;
		    }
		  x1 += mygl.lattice_to_system;
		}
	      z0 += block_size;
	    }
	  y0 += block_size;
	}
      x0 += block_size;
    }
  glEnd ();
  
  
  if (draw_triangles)
    {
      glLineWidth (3.F);
      
      for (boundary_index = 0; boundary_index < BOUNDARIES; boundary_index++)
	{
	  if (boundary_index == INLET_BOUNDARY)
	    {
	      glColor3f (0.5F, 0.5F, 0.5F);
	    }
	  else if (boundary_index == OUTLET_BOUNDARY)
	    {
	      glColor3f (0.0F, 0.5F, 0.0F);
	    }
	  else if (boundary_index == WALL_BOUNDARY)
	    {
	      glColor3f (0.0F, 0.0F, 0.5F);
	    }
	  for (triangle_index = 0; triangle_index < mygl.boundary[ boundary_index ].triangles; triangle_index++)
	    {
	      if (last_triangle.triangle_id == triangle_index &&
		  last_triangle.boundary_id == boundary_index)
		{
		  glBegin (GL_TRIANGLES);
		}
	      else
		{
		  glBegin (GL_LINE_LOOP);
		}
	      triangle_p = &mygl.boundary[ boundary_index ].triangle[ triangle_index ];
	      
	      for (n = 0; n < 3; n++)
		{
		  myglTransformVertex (triangle_p->v[n].pos_x, triangle_p->v[n].pos_y, triangle_p->v[n].pos_z, &x2, &y2, &z2);
		  glVertex3f (x2, y2, z2);
		}
	      glEnd ();
	    }
	}
    }
  
  if (mygl.screen_voxel_coords.i >= 0)
   {
     // vertex quadrant
     pixel_i = mygl.screen_voxel_coords.i;
     pixel_j = mygl.screen_voxel_coords.j;
     
     x1 = (float)pixel_i * (screen.max_x * (2.F / screen_voxels)) - screen.max_x;
     y1 = (float)pixel_j * (screen.max_y * (2.F / screen_voxels)) - screen.max_y;
     
     x2 = x1 + (screen.max_x * (2.F / screen_voxels));
     y2 = y1 + (screen.max_y * (2.F / screen_voxels));
     
     glColor3f (0.F, 0.F, 1.F);
     glBegin (GL_LINE_LOOP);
     
     glVertex3f (x1, y1, 1.0F);
     glVertex3f (x1, y2, 1.0F);
     glVertex3f (x2, y2, 1.0F);
     glVertex3f (x2, y1, 1.0F);
     
     glEnd ();
     
     mygl.screen_voxel_coords.i = -1;
   }
  
  if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
    {
      DisplayTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id);
    }
  else if (passive_mouse_pixel_i >= 0)
    {
      screen_voxel_p = &mygl.screen_to_boundaries_map[ passive_mouse_pixel_i * screen_voxels + passive_mouse_pixel_j ];
      
      DisplaySite (screen_voxel_p->site_i, screen_voxel_p->site_j, screen_voxel_p->site_k);
      
      // site quadrant
      pixel_i = mygl.screen_voxel_coords.i;
      pixel_j = mygl.screen_voxel_coords.j;
      
      x1 = (float)pixel_i * (screen.max_x * (2.F / screen_voxels)) - screen.max_x;
      y1 = (float)pixel_j * (screen.max_y * (2.F / screen_voxels)) - screen.max_y;
      
      x2 = x1 + (screen.max_x * (2.F / screen_voxels));
      y2 = y1 + (screen.max_y * (2.F / screen_voxels));
      
      glColor3f (0.F, 0.F, 1.F);
      glBegin (GL_LINE_LOOP);
      
      glVertex3f (x1, y1, 1.0F);
      glVertex3f (x1, y2, 1.0F);
      glVertex3f (x2, y2, 1.0F);
      glVertex3f (x2, y1, 1.0F);
      
      glEnd ();
    }
  
  glutSwapBuffers ();
}


void DisplaySystemFast (void)
{
  // fast version!
  
  float x1, y1, z1;
  float x2, y2, z2;
  
  int pixel_i, pixel_j;
  int n;
  int triangle_index;
  
  unsigned int boundary_index;
  
  Triangle *triangle_p;
  
  ScreenVoxel *screen_voxel_p;
  
  SuperficialSite *superficial_site_p;
  
  
  display_id = 3;
  
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  
  glColor3f (0.F, 0.F, 0.F);
  
  glLineWidth (3.F);
  glBegin (GL_LINE_LOOP);
  
  myglTransformVertex (-mygl.half_dim_x, -mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (-mygl.half_dim_x, -mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, -mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, -mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  
  glEnd ();
  glBegin (GL_LINE_LOOP);

  myglTransformVertex (-mygl.half_dim_x, +mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (-mygl.half_dim_x, +mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, +mygl.half_dim_y, +mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  myglTransformVertex (+mygl.half_dim_x, +mygl.half_dim_y, -mygl.half_dim_z, &x1, &y1, &z1);	  
  glVertex3f (x1, y1, z1);
  
  glEnd ();
  
  
  for (n = 0; n < screen_voxels * screen_voxels; n++)
    {
      mygl.screen_to_boundaries_map[ n ].site_i = -1;
      mygl.screen_to_boundaries_map[ n ].site_z = 0.F;
    }
  
  glPointSize ((float)point_size);
  glBegin (GL_POINTS);
  
  for (n = 0; n < mygl.superficial_sites; n++)
    {
      superficial_site_p = &mygl.superficial_site[ n ];
      
      x1 = (float)superficial_site_p->i * mygl.lattice_to_system - mygl.half_dim_x;
      y1 = (float)superficial_site_p->j * mygl.lattice_to_system - mygl.half_dim_y;
      z1 = (float)superficial_site_p->k * mygl.lattice_to_system - mygl.half_dim_z;
      
      glColor3f (RedComponent (superficial_site_p->iters), 0.F, 0.F);
      
      myglTransformVertex (x1, y1, z1, &x2, &y2, &z2);
      glVertex3f (x2, y2, z2);
      
      pixel_i = (int)(screen_voxels_screen_max_inv_x * (x2 + screen.max_x));
      pixel_j = (int)(screen_voxels_screen_max_inv_y * (y2 + screen.max_y));
      
      if (pixel_i < 0 || pixel_i >= screen_voxels ||
	  pixel_j < 0 || pixel_j >= screen_voxels)
	{
	  continue;
	}
      screen_voxel_p = &mygl.screen_to_boundaries_map[ pixel_i * screen_voxels + pixel_j ];
      
      if (z2 > screen_voxel_p->site_z)
	{
	  screen_voxel_p->site_z  = z2;
	  screen_voxel_p->site_i = superficial_site_p->i;
	  screen_voxel_p->site_j = superficial_site_p->j;
	  screen_voxel_p->site_k = superficial_site_p->k;
	}
    }
  glEnd ();
  
  
  if (draw_triangles)
    {
      glLineWidth (3.F);
      
      for (boundary_index = 0; boundary_index < BOUNDARIES; boundary_index++)
	{
	  if (boundary_index == INLET_BOUNDARY)
	    {
	      glColor3f (0.5F, 0.5F, 0.5F);
	    }
	  else if (boundary_index == OUTLET_BOUNDARY)
	    {
	      glColor3f (0.0F, 0.5F, 0.0F);
	    }
	  else if (boundary_index == WALL_BOUNDARY)
	    {
	      glColor3f (0.0F, 0.0F, 0.5F);
	    }
	  for (triangle_index = 0; triangle_index < mygl.boundary[ boundary_index ].triangles; triangle_index++)
	    {
	      if (last_triangle.triangle_id == triangle_index &&
		  last_triangle.boundary_id == boundary_index)
		{
		  glBegin (GL_TRIANGLES);
		}
	      else
		{
		  glBegin (GL_LINE_LOOP);
		}
	      triangle_p = &mygl.boundary[ boundary_index ].triangle[ triangle_index ];
	      
	      for (n = 0; n < 3; n++)
		{
		  myglTransformVertex (triangle_p->v[n].pos_x, triangle_p->v[n].pos_y, triangle_p->v[n].pos_z, &x2, &y2, &z2);
		  glVertex3f (x2, y2, z2);
		}
	      glEnd ();
	    }
	}
    }
  

  if (mygl.screen_voxel_coords.i >= 0)
   {
     // vertex quadrant
     pixel_i = mygl.screen_voxel_coords.i;
     pixel_j = mygl.screen_voxel_coords.j;
     
     x1 = (float)pixel_i * (screen.max_x * (2.F / screen_voxels)) - screen.max_x;
     y1 = (float)pixel_j * (screen.max_y * (2.F / screen_voxels)) - screen.max_y;
     
     x2 = x1 + (screen.max_x * (2.F / screen_voxels));
     y2 = y1 + (screen.max_y * (2.F / screen_voxels));
     
     glColor3f (0.F, 0.F, 1.F);
     glBegin (GL_LINE_LOOP);
     
     glVertex3f (x1, y1, 1.0F);
     glVertex3f (x1, y2, 1.0F);
     glVertex3f (x2, y2, 1.0F);
     glVertex3f (x2, y1, 1.0F);
     
     glEnd ();
     
     mygl.screen_voxel_coords.i = -1;
   }
  
  if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
    {
      DisplayTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id);
    }
  else if (passive_mouse_pixel_i >= 0)
    {
      screen_voxel_p = &mygl.screen_to_boundaries_map[ passive_mouse_pixel_i * screen_voxels + passive_mouse_pixel_j ];
      
      DisplaySite (screen_voxel_p->site_i, screen_voxel_p->site_j, screen_voxel_p->site_k);
      
      // site quadrant
      pixel_i = mygl.screen_voxel_coords.i;
      pixel_j = mygl.screen_voxel_coords.j;
      
      x1 = (float)pixel_i * (screen.max_x * (2.F / screen_voxels)) - screen.max_x;
      y1 = (float)pixel_j * (screen.max_y * (2.F / screen_voxels)) - screen.max_y;
      
      x2 = x1 + (screen.max_x * (2.F / screen_voxels));
      y2 = y1 + (screen.max_y * (2.F / screen_voxels));
      
      glColor3f (0.F, 0.F, 1.F);
      glBegin (GL_LINE_LOOP);
      
      glVertex3f (x1, y1, 1.0F);
      glVertex3f (x1, y2, 1.0F);
      glVertex3f (x2, y2, 1.0F);
      glVertex3f (x2, y1, 1.0F);
      
      glEnd ();
    }
  
  glutSwapBuffers ();
}


void GLUTCALLBACK Display (void)
{
  void (*DisplayPointer[4])(void);
  
  DisplayPointer[0] = DisplayNull;
  DisplayPointer[1] = DisplaySlice;
  DisplayPointer[2] = DisplaySystemFast;
  DisplayPointer[3] = DisplaySystemSlow;
  
  (*DisplayPointer[ display_id ]) ();
}


void ReadFirstSlice (int *pixels_x, int *pixels_y, int *voxel_size)
{
  FILE *temp_file_ptr, *image_file_ptr;
  
  char partial_file_name[256], image_file_name_input[256], *image_file_name_output;
  char *partial_command, total_command[256];
  
  
  temp_file_ptr = fopen ("temp_file.txt", "r");
  
  fscanf (temp_file_ptr, "%s\n", partial_file_name);
  fclose (temp_file_ptr);
  
  sprintf (image_file_name_input, "%s%s", input_path, partial_file_name);
  
  partial_command = "convert";
  image_file_name_output = "temp_file.ppm";
  
  sprintf (total_command, "%s %s %s",
	   partial_command, image_file_name_input, image_file_name_output);
  system (total_command);
  
  image_file_ptr = fopen (image_file_name_output, "r");
  
  fscanf (image_file_ptr, "%s\n%i %i\n%i\n", ppm_type, pixels_x, pixels_y, voxel_size);
  
  printf ("input image pixels x, y: %i %i, voxel size: %i\n",
	  *pixels_x, *pixels_y, *voxel_size);
  fflush (stdout);
  
  fclose (image_file_ptr);
}


void ReadSlice255 (int slice_id, int gray_max, unsigned char slice_row_data[])
{
  FILE *temp_file_ptr, *image_file_ptr;
  
  char partial_file_name[256], image_file_name_input[256], *image_file_name_output;
  char *partial_command, total_command[256];
  
  int pixels_x, pixels_y;
  int i, j;
  
  
  printf ("processing slice %i ...\n", slice_id);
  fflush (stdout);
  
  temp_file_ptr = fopen ("temp_file.txt", "r");
  
  for (i = 0; i < slice_id; i++)
    {
      fscanf (temp_file_ptr, "%s\n", partial_file_name);
    }
  fclose (temp_file_ptr);
  
  sprintf (image_file_name_input, "%s%s", input_path, partial_file_name);
  
  partial_command = "convert";
  image_file_name_output = "temp_file.ppm";
  
  sprintf (total_command, "%s %s %s",
	   partial_command, image_file_name_input, image_file_name_output);
  system (total_command);
  
  
  if (strcmp (ppm_type, "P3"))
    {
      image_file_ptr = fopen (image_file_name_output, "rb");
      
      fscanf (image_file_ptr, "%s\n%i %i\n%i\n", ppm_type, &pixels_x, &pixels_y, &gray_max);
      
      for (j = 0; j < mygl.input_image_pix_y; j++)
	{
	  fread (slice_row_data, sizeof(unsigned char), mygl.input_image_pix_x * 3, image_file_ptr);
	  
	  for (i = 0; i < mygl.input_image_pix_x; i++)
	    {
	      mygl.medical_data[ (i*mygl.input_image_pix_y+j)*mygl.input_slices+slice_id ] =
		slice_row_data[ i*3 ];
	    }
	}
      fclose (image_file_ptr);
    }
  else
    {
      image_file_ptr = fopen (image_file_name_output, "r");
      
      ifstream inputFile;
      inputFile.open(image_file_name_output);
      
      inputFile >> ppm_type;
      inputFile >> pixels_x;
      inputFile >> pixels_y;
      inputFile >> gray_max;
      
      for (j = 0; j < mygl.input_image_pix_y; j++)
	{
	  for (i = 0; i < mygl.input_image_pix_x * 3; i++)
	    {
	      int temp;
	      
	      inputFile >> temp;
	      slice_row_data[i] = temp;
	    }
	  for (i = 0; i < mygl.input_image_pix_x; i++)
	    {
	      mygl.medical_data[ (i*mygl.input_image_pix_y+j)*mygl.input_slices+slice_id ] =
		slice_row_data[ i*3 ];
	    }
	}
      inputFile.close();
    }
}


void ReadSlice65536 (int slice_id, int gray_max, unsigned short int slice_row_data[])
{
  FILE *temp_file_ptr, *image_file_ptr;
  
  char partial_file_name[256], image_file_name_input[256], *image_file_name_output;
  char *partial_command, total_command[256];
  
  int pixels_x, pixels_y;
  int i, j;
  
  
  printf ("processing slice %i ...\n", slice_id);
  fflush (stdout);
  
  temp_file_ptr = fopen ("temp_file.txt", "r");
  
  for (i = 0; i < slice_id; i++)
    {
      fscanf (temp_file_ptr, "%s\n", partial_file_name);
    }
  fclose (temp_file_ptr);
  
  sprintf (image_file_name_input, "%s%s", input_path, partial_file_name);
  
  partial_command = "convert";
  image_file_name_output = "temp_file.ppm";
  
  sprintf (total_command, "%s %s %s",
	   partial_command, image_file_name_input, image_file_name_output);
  system (total_command);
  
  
  if (strcmp (ppm_type, "P3"))
    {
      image_file_ptr = fopen (image_file_name_output, "rb");
      
      fscanf (image_file_ptr, "%s\n%i %i\n%i\n", ppm_type, &pixels_x, &pixels_y, &gray_max);
      
      for (j = 0; j < mygl.input_image_pix_y; j++)
	{
	  fread (slice_row_data, sizeof(unsigned short int), mygl.input_image_pix_x * 3, image_file_ptr);
	  
	  for (i = 0; i < mygl.input_image_pix_x; i++)
	    {
	      mygl.medical_data[ (i*mygl.input_image_pix_y+j)*mygl.input_slices+slice_id ] =
		slice_row_data[ i*3 ];
	    }
	}
      fclose (image_file_ptr);
    }
  else
    {
      image_file_ptr = fopen (image_file_name_output, "r");
      
      ifstream inputFile;
      inputFile.open(image_file_name_output);
      
      inputFile >> ppm_type;
      inputFile >> pixels_x;
      inputFile >> pixels_y;
      inputFile >> gray_max;
      
      for (j = 0; j < mygl.input_image_pix_y; j++)
	{
	  for (i = 0; i < mygl.input_image_pix_x * 3; i++)
	    {
	      int temp;
	      
	      inputFile >> temp;
	      slice_row_data[i] = temp;
	    }
	  for (i = 0; i < mygl.input_image_pix_x; i++)
	    {
	      mygl.medical_data[ (i*mygl.input_image_pix_y+j)*mygl.input_slices+slice_id ] =
		slice_row_data[ i*3 ];
	    }
	}
      inputFile.close();
    }
}


void ReadConfig (void)
{
  FILE *temp_file_ptr;
  
  float seconds;
  
  int pixels_x, pixels_y, gray_max;
  int i;
  
  char *first_command_part, *last_command_part;
  char my_command[256];
  
  
  // the names of the files corresponding to the data slices are
  // written in a temporary file
  
  first_command_part = "ls";
  last_command_part  = "| less | wc > temp_file.txt";
  
  sprintf (my_command, "%s %s %s", first_command_part, input_path, last_command_part);
  system (my_command);
  
  temp_file_ptr = fopen ("temp_file.txt", "r");
  fscanf (temp_file_ptr, "%i ", &mygl.input_slices);
  fclose (temp_file_ptr);
  
  temp_file_ptr = fopen ("temp_file.txt", "r");
  
  first_command_part = "ls -1";
  last_command_part  = "> temp_file.txt";
  
  sprintf (my_command, "%s %s %s", first_command_part, input_path, last_command_part);
  system (my_command);
  
  fclose (temp_file_ptr);
  
  
  ReadFirstSlice (&pixels_x, &pixels_y, &gray_max);
  
  mygl.input_image_pix_x = pixels_x;
  mygl.input_image_pix_y = pixels_y;
  
  mygl.output_image_pix_x = (int)(mygl.input_image_pix_x * res_factor);
  mygl.output_image_pix_y = (int)(mygl.input_image_pix_y * res_factor);
  
  mygl.output_slices = (int)(mygl.input_slices * (slice_size / pixel_size) *
			     (float)mygl.output_image_pix_x / (float)mygl.input_image_pix_x);
  
  mygl.scale_x = res_factor;
  mygl.scale_y = res_factor;
  mygl.scale_z = res_factor * slice_size / pixel_size;
  
  mygl.scale_inv_x = 1.F / mygl.scale_x;
  mygl.scale_inv_y = 1.F / mygl.scale_y;
  mygl.scale_inv_z = 1.F / mygl.scale_z;
  
  printf (" output_slices: %i\n", mygl.output_slices);
  fflush (stdout);
  
  // system parameters setup
  
  mygl.lattice_to_system = 1.F;
  
  mygl.block_size = 8;
  mygl.shift = 3;
  
  mygl.blocks_x = mygl.output_image_pix_x >> mygl.shift;
  mygl.blocks_y = mygl.output_image_pix_y >> mygl.shift;
  mygl.blocks_z = mygl.output_slices      >> mygl.shift;
  
  if ((mygl.blocks_x << mygl.shift) < mygl.output_image_pix_x) ++mygl.blocks_x;
  if ((mygl.blocks_y << mygl.shift) < mygl.output_image_pix_y) ++mygl.blocks_y;
  if ((mygl.blocks_z << mygl.shift) < mygl.output_slices     ) ++mygl.blocks_z;
  
  mygl.sites_x = mygl.blocks_x * mygl.block_size;
  mygl.sites_y = mygl.blocks_y * mygl.block_size;
  mygl.sites_z = mygl.blocks_z * mygl.block_size;
  
  mygl.dim_x = (float)mygl.sites_x * mygl.lattice_to_system;
  mygl.dim_y = (float)mygl.sites_y * mygl.lattice_to_system;
  mygl.dim_z = (float)mygl.sites_z * mygl.lattice_to_system;
  
  mygl.half_dim_x = 0.5F * mygl.dim_x;
  mygl.half_dim_y = 0.5F * mygl.dim_y;
  mygl.half_dim_z = 0.5F * mygl.dim_z;
  
  mygl.system_size = fmaxf(mygl.dim_x, fmaxf(mygl.dim_y, mygl.dim_z));
  
  mygl.sites_in_a_block = mygl.block_size * mygl.block_size * mygl.block_size;
  
  mygl.blocks = mygl.blocks_x * mygl.blocks_y * mygl.blocks_z;
  
  
  // initial setup of the medical and system datasets
  
  mygl.medical_data = (unsigned short int *)malloc(sizeof(unsigned short int) * mygl.input_slices * pixels_x * pixels_y);
  
  mygl.data_block = (DataBlock *)malloc(sizeof(DataBlock) * mygl.blocks);
  
  for (i = 0; i < mygl.blocks; i++)
    {
      mygl.data_block[ i ].site_data = NULL;
      mygl.data_block[ i ].site_label = NULL;
      mygl.data_block[ i ].site_iters = NULL;
    }
  
  if (gray_max == 255)
    {
      mygl.slice_row_data = (unsigned char *)malloc(sizeof(unsigned short int) * pixels_x * 3);
    }
  else
    {
      mygl.slice_row_data = (unsigned short int *)malloc(sizeof(unsigned char) * pixels_x * 3);
    }
  
  // the medical dataset is filtered and stored in a bi-level data structure
  
  seconds = myClock ();
  
  for (i = 0; i < mygl.input_slices; i++)
    {
      if (gray_max == 255)
	{
	  ReadSlice255 (i, gray_max, (unsigned char *)mygl.slice_row_data);
	}
      else
	{
	  ReadSlice65536 (i, gray_max, (unsigned short int *)mygl.slice_row_data);
	}
    }
  printf ("seconds to read: %.3f\n", myClock () - seconds);
  fflush (stdout);
  
  free(mygl.slice_row_data);
  
  system ("rm temp_file.txt");
  system ("rm temp_file.ppm");
  
  
  mygl.superficial_sites_max = 100000;
  mygl.superficial_site = (SuperficialSite *)malloc(sizeof(SuperficialSite) * mygl.superficial_sites_max);
  
  mygl.stored_blocks_max = 1000;
  mygl.stored_block = (int *)malloc(sizeof(int) * mygl.stored_blocks_max);
  
  site_location_a = (SiteLocation *)malloc(sizeof(SiteLocation) * BUFFERS_SIZE);
  site_location_b = (SiteLocation *)malloc(sizeof(SiteLocation) * BUFFERS_SIZE);
}


void WriteCheckpoint (char *file_name)
{
  FILE *system_config;
  XDR xdr_config;
  
  double lattice_to_system;
  
  int i, m, n;
  
  unsigned int data;
  
  
  system_config = fopen (file_name, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
  
  xdr_float (&xdr_config, &slice_size);
  xdr_float (&xdr_config, &pixel_size);
  xdr_float (&xdr_config, &res_factor);
  
  xdr_int     (&xdr_config, &mygl.input_image_pix_x);
  xdr_int     (&xdr_config, &mygl.input_image_pix_y);
  xdr_int     (&xdr_config, &mygl.input_slices);
  
  xdr_int     (&xdr_config, &mygl.selected_pixel_x);
  xdr_int     (&xdr_config, &mygl.selected_pixel_y);
  xdr_int     (&xdr_config, &mygl.selected_slice);
  xdr_u_short (&xdr_config, &mygl.selected_gray);
  
  lattice_to_system = (double)mygl.lattice_to_system;
  
  xdr_double (&xdr_config, &lattice_to_system);
  xdr_int    (&xdr_config, &mygl.blocks_x);
  xdr_int    (&xdr_config, &mygl.blocks_y);
  xdr_int    (&xdr_config, &mygl.blocks_z);
  xdr_int    (&xdr_config, &mygl.block_size);
  

  for (n = 0; n < mygl.input_slices * mygl.input_image_pix_x * mygl.input_image_pix_y; n++)
    {
      data = mygl.medical_data[ n ];
      
      xdr_u_int (&xdr_config, &data);
    }
  
  for (n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &mygl.boundary[ n ].triangles);
      
      for (m = 0; m < mygl.boundary[ n ].triangles; m++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].v[ i ].pos_x);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].v[ i ].pos_y);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].v[ i ].pos_z);
	    }
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].pressure_avg);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].pressure_amp);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].pressure_phs);
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void ReadCheckpoint (char *file_name)
{
  FILE *system_config;
  XDR xdr_config;
  
  double lattice_to_system;
  
  int pixels_x, pixels_y;
  int i, m, n;
  
  unsigned int data;
  
  
  system_config = fopen (file_name, "r");
  xdrstdio_create (&xdr_config, system_config, XDR_DECODE);
  
  xdr_float (&xdr_config, &slice_size);
  xdr_float (&xdr_config, &pixel_size);
  xdr_float (&xdr_config, &res_factor);
  
  xdr_int     (&xdr_config, &mygl.input_image_pix_x);
  xdr_int     (&xdr_config, &mygl.input_image_pix_y);
  xdr_int     (&xdr_config, &mygl.input_slices);
  
  xdr_int     (&xdr_config, &mygl.selected_pixel_x);
  xdr_int     (&xdr_config, &mygl.selected_pixel_y);
  xdr_int     (&xdr_config, &mygl.selected_slice);
  xdr_u_short (&xdr_config, &mygl.selected_gray);
  
  xdr_double (&xdr_config, &lattice_to_system);
  xdr_int    (&xdr_config, &mygl.blocks_x);
  xdr_int    (&xdr_config, &mygl.blocks_y);
  xdr_int    (&xdr_config, &mygl.blocks_z);
  xdr_int    (&xdr_config, &mygl.block_size);
  
  pixels_x = mygl.input_image_pix_x;
  pixels_y = mygl.input_image_pix_y;
  
  mygl.output_image_pix_x = (int)(mygl.input_image_pix_x * res_factor);
  mygl.output_image_pix_y = (int)(mygl.input_image_pix_y * res_factor);
  
  mygl.output_slices = (int)(mygl.input_slices * (slice_size / pixel_size) *
			     (float)mygl.output_image_pix_x / (float)mygl.input_image_pix_x);
  
  mygl.scale_x = (float)mygl.output_image_pix_x / (float)mygl.input_image_pix_x;
  mygl.scale_y = (float)mygl.output_image_pix_y / (float)mygl.input_image_pix_y;
  mygl.scale_z = (float)mygl.output_slices / (float)mygl.input_slices;
  
  mygl.scale_inv_x = 1.F / mygl.scale_x;
  mygl.scale_inv_y = 1.F / mygl.scale_y;
  mygl.scale_inv_z = 1.F / mygl.scale_z;
  
  mygl.lattice_to_system = lattice_to_system;
  
  mygl.sites_x = mygl.blocks_x * mygl.block_size;
  mygl.sites_y = mygl.blocks_y * mygl.block_size;
  mygl.sites_z = mygl.blocks_z * mygl.block_size;
  
  mygl.dim_x = (float)mygl.sites_x * lattice_to_system;
  mygl.dim_y = (float)mygl.sites_y * lattice_to_system;
  mygl.dim_z = (float)mygl.sites_z * lattice_to_system;
  
  mygl.half_dim_x = 0.5F * mygl.dim_x;
  mygl.half_dim_y = 0.5F * mygl.dim_y;
  mygl.half_dim_z = 0.5F * mygl.dim_z;
  
  mygl.system_size = fmaxf(mygl.dim_x, fmaxf(mygl.dim_y, mygl.dim_z));
  
  mygl.sites_in_a_block = mygl.block_size * mygl.block_size * mygl.block_size;
  
  mygl.blocks = mygl.blocks_x * mygl.blocks_y * mygl.blocks_z;
  
  i = mygl.block_size;
  
  mygl.shift = 0;

  while (i > 1)
    {
      i >>= 1;

      ++mygl.shift;
    }
  
  // initial setup of the medical and system datasets
  
  mygl.medical_data = (unsigned short int *)malloc(sizeof(unsigned short int) * mygl.input_slices * pixels_x * pixels_y);
  
  for (n = 0; n < mygl.input_slices * pixels_x * pixels_y; n++)
    {
      xdr_u_int (&xdr_config, &data);
      
      mygl.medical_data[ n ] = data;
    }
  
  mygl.data_block = (DataBlock *)malloc(sizeof(DataBlock) * mygl.blocks);
  
  for (n = 0; n < mygl.blocks; n++)
    {
      mygl.data_block[ n ].site_data = NULL;
      mygl.data_block[ n ].site_iters = NULL;
      mygl.data_block[ n ].is_void = 1;
    }
  
  
  mygl.boundary[ INLET_BOUNDARY  ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  mygl.boundary[ OUTLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  mygl.boundary[ WALL_BOUNDARY   ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  
  for (n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &mygl.boundary[ n ].triangles);
      
      for (m = 0; m < mygl.boundary[ n ].triangles; m++)
	{
	  for (i = 0; i < 3; i++)
	    {
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].v[ i ].pos_x);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].v[ i ].pos_y);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].v[ i ].pos_z);
	    }
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].pressure_avg);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].pressure_amp);
	      xdr_float (&xdr_config, &mygl.boundary[ n ].triangle[ m ].pressure_phs);
	    }
	}
    }
  xdr_destroy (&xdr_config);
  
  
  mygl.screen_to_boundaries_map = (ScreenVoxel *)malloc(sizeof(ScreenVoxel) * screen_voxels * screen_voxels);
  
  mygl.superficial_sites_max = 100000;
  mygl.superficial_site = (SuperficialSite *)malloc(sizeof(SuperficialSite) * mygl.superficial_sites_max);
  
  mygl.stored_blocks_max = 1000;
  mygl.stored_block = (int *)malloc(sizeof(int) * mygl.stored_blocks_max);
  
  site_location_a = (SiteLocation *)malloc(sizeof(SiteLocation) * BUFFERS_SIZE);
  site_location_b = (SiteLocation *)malloc(sizeof(SiteLocation) * BUFFERS_SIZE);
  
  
  myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
  
  myglFluidSitesIterativeSearching (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
  
  myglSetBoundaryConfigurations ();
}


void WriteConfig (char *file_name)
{
  FILE *system_config;
  XDR xdr_config;

  double lattice_to_system;
  
  int i, j, k;
  int m, n;
  int are_all_solid_sites;
  int flag;
  
  unsigned int site_data;

  DataBlock *data_block_p;
  
  
  system_config = fopen (file_name, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
  
  lattice_to_system = (double)mygl.lattice_to_system;
  
  xdr_double (&xdr_config, &lattice_to_system);
  xdr_int    (&xdr_config, &mygl.blocks_x);
  xdr_int    (&xdr_config, &mygl.blocks_y);
  xdr_int    (&xdr_config, &mygl.blocks_z);
  xdr_int    (&xdr_config, &mygl.block_size);
  
  n = -1;
  
  for (i = 0; i < mygl.blocks_x; i++)
    {
      for (j = 0; j < mygl.blocks_y; j++)
	{
	  for (k = 0; k < mygl.blocks_z; k++)
	    {
	      data_block_p = &mygl.data_block[ ++n ];
	      
	      flag = 0;
	      
	      if (data_block_p->site_data == NULL)
		{
		  xdr_int (&xdr_config, &flag);
		  continue;
		}
	      are_all_solid_sites = 1;
	      
	      for (m = 0; m < mygl.sites_in_a_block; m++)
		{
		  if (data_block_p->site_data[ m ] != SOLID_TYPE &&
		      data_block_p->site_data[ m ] != NULL_TYPE)
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
	      
	      for (m = 0; m < mygl.sites_in_a_block; m++)
		{
		  site_data = data_block_p->site_data[ m ];
		  
		  if (site_data == NULL_TYPE)
		    {
		      site_data = SOLID_TYPE;
		    }
		  xdr_u_int (&xdr_config, &site_data);
		}
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void WritePars (char *file_name)
{
  FILE *pars = fopen (file_name, "w");
  
  float nx, ny, nz;
  
  int n;
  
  Triangle *triangle_p;
  
  
  fprintf (pars, "%i\n", mygl.boundary[ INLET_BOUNDARY ].triangles);
  
  for (n = 0; n < mygl.boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      fprintf (pars, "%f %f %f\n",
	       mygl.boundary[ INLET_BOUNDARY ].triangle[n].pressure_avg,
	       mygl.boundary[ INLET_BOUNDARY ].triangle[n].pressure_amp,
	       mygl.boundary[ INLET_BOUNDARY ].triangle[n].pressure_phs);
    }
  
  fprintf (pars, "%i\n", mygl.boundary[ OUTLET_BOUNDARY ].triangles);
  
  for (n = 0; n < mygl.boundary[ OUTLET_BOUNDARY ].triangles; n++)
    {
      fprintf (pars, "%f %f %f\n",
	       mygl.boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_avg,
	       mygl.boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_amp,
	       mygl.boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_phs);
    }
  for (n = 0; n < mygl.boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      triangle_p = &mygl.boundary[ INLET_BOUNDARY ].triangle[n];
      
      myglTriangleNormal (triangle_p->v[0].pos_x, triangle_p->v[0].pos_y, triangle_p->v[0].pos_z,
			  triangle_p->v[1].pos_x, triangle_p->v[1].pos_y, triangle_p->v[1].pos_z,
			  triangle_p->v[2].pos_x, triangle_p->v[2].pos_y, triangle_p->v[2].pos_z,
			  &nx, &ny, &nz);
      
      fprintf (pars, "%f %f %f\n", nx, ny, nz);
    }
  fclose (pars);
}


void myglEnd (void)
{
  free(site_location_a);
  free(site_location_b);
  
  
  free(mygl.stored_block);
  mygl.stored_block = NULL;
  
  
  free(mygl.superficial_site);
  mygl.superficial_site = NULL;
  
  
  for (int i = 0; i < mygl.blocks; i++)
    {
      if (mygl.data_block[ i ].site_data != NULL)
	{
	  free(mygl.data_block[ i ].site_data);
	  mygl.data_block[ i ].site_data = NULL;
	}
    }
  free(mygl.data_block);
  mygl.data_block = NULL;
  
  
  free(mygl.medical_data);
}


void myglInitBoundaries (void)
{
  mygl.boundary[ INLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  mygl.boundary[ INLET_BOUNDARY ].triangles = 0;
  
  mygl.boundary[ OUTLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  mygl.boundary[ OUTLET_BOUNDARY ].triangles = 0;
  
  mygl.boundary[ WALL_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U << BOUNDARY_ID_BITS));
  mygl.boundary[ WALL_BOUNDARY ].triangles = 0;
  
  
  mygl.screen_to_boundaries_map = (ScreenVoxel *)malloc(sizeof(ScreenVoxel) * screen_voxels * screen_voxels);
}

void myglEndBoundaries (void)
{
  free(mygl.screen_to_boundaries_map);
  mygl.screen_to_boundaries_map = NULL;

  
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      free(mygl.boundary[ n ].triangle);
      mygl.boundary[ n ].triangle = NULL;
      
      mygl.boundary[ n ].triangles = 0;
    }
}


void  MoveSceneCenter (float t_x, float t_y)
{
  float x, y, z;
  
  
  myglTransformVertex (scene_center_x, scene_center_y, scene_center_z, &x, &y, &z);
      
  x += t_x;
  y += t_y;
  
  myglAntiTransformVertex (x, y, z, &scene_center_x, &scene_center_y, &scene_center_z);
}


void GLUTCALLBACK KeybordFunction (unsigned char key, int x, int y)
{
  int triangle_id;
  
  ScreenVoxel *screen_voxel_p;
  
  
  if (key == 'S')
    {
      printf("Opening ppm file: ./image.ppm\n");
      SaveWindowImage ("./image.ppm");
      
      if (display_id != 1)
	{
	  printf("Opening output pars file: %s\n", output_pars_name);
	  WritePars (output_pars_name);
	  
	  printf("Opening output config file: %s\n", output_config_name);
	  WriteConfig (output_config_name);
	}
      printf("Opening checkpoint file: %s\n", checkpoint_name);
      WriteCheckpoint (checkpoint_name);
    }
  else if (key == 'q')
    {
      printf("Opening ppm file: ./image.ppm\n");
      SaveWindowImage ("./image.ppm");
      
      if (display_id != 1)
	{
	  printf("Opening output pars file: %s\n", output_pars_name);
	  WritePars (output_pars_name);
	  
	  printf("Opening output config file: %s\n", output_config_name);
	  WriteConfig (output_config_name);
	}
      printf("Opening checkpoint file: %s\n", checkpoint_name);
      WriteCheckpoint (checkpoint_name);
      
      myglEndBoundaries ();
      
      myglEnd ();
      
      exit(0);
    }
  else if (key == 'm')
    {
      display_id = 1;
    }
  
  if (display_id == 1)
    {
      if (key == '>')
	{
	  ++mygl.selected_slice;
	  mygl.selected_slice = min(mygl.input_slices - 1, mygl.selected_slice);
	  printf (" slice selected: %i\n", mygl.selected_slice);
	}
      else if (key == '<')
	{
	  --mygl.selected_slice;
	  mygl.selected_slice = max(0, mygl.selected_slice);
	  printf (" slice selected: %i\n", mygl.selected_slice);
	}
      return;
    }
  
  if (key == 'c')
    {
      ortho_x = 0.5F * mygl.system_size;
      ortho_y = 0.5F * mygl.system_size;
      
      longitude = 0.F;
      latitude = 0.F;
      viewpoint_radius = 10.F * mygl.system_size;
      
      zoom = 1.5F;
      
      scene_center_x = 0.F;
      scene_center_y = 0.F;
      scene_center_z = 0.F;
      
      Projection ();
    }
  else if (key == 'z')
    {
      if (last_triangle.triangle_id == -1)
	{
	  zoom *= 1.1F;
	  Projection ();
	}
      else
	{
	  EditLastTriangle (0.F, 0.F, 0.F, 0.F, 1.1F);
	}
    }
  else if (key == 'Z')
    {
      if (last_triangle.triangle_id == -1)
	{
	  zoom /= 1.1F;
	  Projection ();
	}
      else
	{
	  EditLastTriangle (0.F, 0.F, 0.F, 0.F, 1.F / 1.1F);
	}
    }
  else if (key == 'a')
    {
      if (last_triangle.triangle_id == -1)
	{
	  MoveSceneCenter (5.F * mygl.lattice_to_system, 0.F);
	  Projection ();
	}
      else
	{
	  EditLastTriangle (-0.5F * mygl.lattice_to_system, 0.F, 0.F, 0.F, 1.0F);
	}
    }
  else if (key == 'd')
    {
      if (last_triangle.triangle_id == -1)
	{
	  MoveSceneCenter (-5.F * mygl.lattice_to_system, 0.F);
	  Projection ();
	}
      else
	{
	  EditLastTriangle (0.5F * mygl.lattice_to_system, 0.F, 0.F, 0.F, 1.0F);
	}
    }
  else if (key == 's')
    {
      if (last_triangle.triangle_id == -1)
	{
	  MoveSceneCenter (0.F, 5.F * mygl.lattice_to_system);
	  Projection ();
	}
      else
	{
	  EditLastTriangle (0.F, -0.5F * mygl.lattice_to_system, 0.F, 0.F, 1.0F);
	}
    }
  else if (key == 'w')
    {
      if (last_triangle.triangle_id == -1)
	{
	  MoveSceneCenter (0.F, -5.F * mygl.lattice_to_system);
	  Projection ();
	}
      else
	{
	  EditLastTriangle (0.F, 0.5F * mygl.lattice_to_system, 0.F, 0.F, 1.0F);
	}
    }
  else if (key == 'D')
    {
      if (last_triangle.triangle_id != -1)
	{
	  DeleteLastTriangle ();
	  
	  last_triangle.triangle_id = -1;
	  last_triangle.boundary_id = -1;
	}
    }
  else if (key == '1')
    {
      if (passive_mouse_pixel_i == -1)
	{
	  return;
	}
      screen_voxel_p = &mygl.screen_to_boundaries_map[ passive_mouse_pixel_i * screen_voxels +
						       passive_mouse_pixel_j ];
      
      triangle_id = myglCreateOptimizedTriangle (INLET_BOUNDARY,
						 screen_voxel_p->site_i, screen_voxel_p->site_j, screen_voxel_p->site_k);
      
      if (triangle_id != -1)
	{
	  last_triangle.triangle_id = triangle_id;
	  last_triangle.boundary_id = INLET_BOUNDARY;
	}
    }
  else if (key == '2')
    {
     if (passive_mouse_pixel_i == -1)
	{
	  return;
	}
      screen_voxel_p = &mygl.screen_to_boundaries_map[ passive_mouse_pixel_i * screen_voxels +
						       passive_mouse_pixel_j ];
      
      triangle_id = myglCreateOptimizedTriangle (OUTLET_BOUNDARY,
						 screen_voxel_p->site_i, screen_voxel_p->site_j, screen_voxel_p->site_k);
      
      if (triangle_id != -1)
	{
	  last_triangle.triangle_id = triangle_id;
	  last_triangle.boundary_id = OUTLET_BOUNDARY;
	}
    }
  else if (key == '3')
    {
     if (passive_mouse_pixel_i == -1)
	{
	  return;
	}
      screen_voxel_p = &mygl.screen_to_boundaries_map[ passive_mouse_pixel_i * screen_voxels +
						       passive_mouse_pixel_j ];
      
      triangle_id = myglCreateOptimizedTriangle (WALL_BOUNDARY,
						 screen_voxel_p->site_i, screen_voxel_p->site_j, screen_voxel_p->site_k);
    
      if (triangle_id != -1)
	{
	  last_triangle.triangle_id = triangle_id;
	  last_triangle.boundary_id = WALL_BOUNDARY;
	}
    }
  else if (key == 't')
    {
      draw_triangles = !draw_triangles;
    }
  else if (key == 'p')
    {
      point_size = min(5, ++point_size);
    }
  else if (key == 'P')
    {
      point_size = max(1, --point_size);
    }
  if (key == 'g')
    {
      mygl.selected_gray += 1;
      mygl.selected_gray = min(1 << 16, mygl.selected_gray);
      
      myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
    }
  else if (key == 'G')
    {
      mygl.selected_gray -= 1;
      mygl.selected_gray = max(0, mygl.selected_gray);
      
      myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
    }
  if (key == 'h')
    {
      mygl.selected_gray += 1 << 4;
      mygl.selected_gray = min(1 << 16, mygl.selected_gray);
      
      myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
    }
  else if (key == 'H')
    {
      mygl.selected_gray -= 1 << 4;
      mygl.selected_gray = max(0, mygl.selected_gray);
      
      myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
    }
  else if (key == 'j')
    {
      if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
	{
	  ChangeTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id, 1.0, 0., 0.);
	}
    }
  else if (key == 'J')
    {
      if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
	{
	  ChangeTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id, -1.0, 0., 0.);
	}
    }
  else if (key == 'k')
    {
      if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
	{
	  ChangeTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id, 0., 0.1, 0.);
	}
    }
  else if (key == 'K')
    {
      if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
	{
	  ChangeTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id, 0., -0.1, 0.);
	}
    }
  else if (key == 'l')
    {
      if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
	{
	  ChangeTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id, 0., 0., 1.);
	}
    }
  else if (key == 'L')
    {
      if (last_triangle.boundary_id == INLET_BOUNDARY || last_triangle.boundary_id == OUTLET_BOUNDARY)
	{
	  ChangeTrianglePars (last_triangle.boundary_id, last_triangle.triangle_id, 0., 0., -1.);
	}
    }
  else if (key == 'r')
    {
      res_factor += 1.;
      printf ("res factor: %.1f\n", res_factor);
      
      myglRescaleSystemResolution ();
      myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
      myglRescaleTriangles (res_factor / (res_factor - 1.F));
      
      myglRescaleViewpoint (res_factor / (res_factor - 1.F));
    }
  else if (key == 'R')
    {
      if (res_factor >= 1.999F)
	{
	  res_factor -= 1.F;
	  printf ("res factor: %.1f\n", res_factor);
	  
	  myglRescaleSystemResolution ();
	  myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
	  myglRescaleTriangles (res_factor / (res_factor + 1.F));
	  
	  myglRescaleViewpoint (res_factor / (res_factor + 1.F));
	}
    }
  display_id = 2;
}

void GLUTCALLBACK SpecialKeybordFunction (int key, int x, int y)
{
  float longitude_inc = 0.F;
  float latitude_inc  = 0.F;
  
  
  if (display_id == 1)
    {
      return;
    }
  
  if (key ==  GLUT_KEY_LEFT)
    {
      longitude_inc = -2.F;
    }
  else if (key == GLUT_KEY_RIGHT)
    {
      longitude_inc = 2.F;
    }
  else if (key == GLUT_KEY_DOWN)
    {
      latitude_inc = -2.F;
    }
  else if (key == GLUT_KEY_UP)
    {
      latitude_inc = 2.F;
    }
  
  if (last_triangle.triangle_id == -1)
    {
      longitude -= longitude_inc;
      latitude  -= latitude_inc;
      
      Projection ();
    }
  else
    {
      EditLastTriangle (0.F, 0.F, longitude_inc * DEG_TO_RAD, latitude_inc * DEG_TO_RAD, 1.F);
    }
  display_id = 2;
}

void GLUTCALLBACK MouseFunction (int button, int state, int x0, int y0)
{
  float x1, y1, z1;
  float x2, y2, z2;
  
  int mouse_pixel_i, mouse_pixel_j;
  int vertex_pixel_i, vertex_pixel_j;
  int i, j, k, n;
  int triangle_index;
  
  unsigned int boundary_index;
  
  Triangle *triangle_p;
  
  ScreenVoxel *screen_voxel_p;
  
  
  y0 = viewport_pixels_y - y0;
  
  passive_mouse_pixel_i = -1;
  passive_mouse_pixel_j = -1;
  
  if (display_id == 1)
    {
      if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
	  mygl.selected_pixel_x = max(0, min(mygl.input_image_pix_x - 1, (int)((float)(x0 * mygl.input_image_pix_x) / viewport_pixels_x)));
	  mygl.selected_pixel_y = max(0, min(mygl.input_image_pix_y - 1, (int)((float)(y0 * mygl.input_image_pix_y) / viewport_pixels_y)));
	  
	  myglReconstructSystem (mygl.selected_pixel_x, mygl.selected_pixel_y, mygl.selected_slice, mygl.selected_gray);
	  
	  display_id = 2;
	}
      return;
    }
  //display_id = 0;
  
  mouse_pixel_i = (int)(screen_voxels * (float)x0 / (float)viewport_pixels_x);
  mouse_pixel_j = (int)(screen_voxels * (float)y0 / (float)viewport_pixels_y);
  
  if (mouse_pixel_i < 0 || mouse_pixel_i >= screen_voxels ||
      mouse_pixel_j < 0 || mouse_pixel_j >= screen_voxels)
    {
      mygl.screen_voxel_pointed = NULL;
      return;
    }
  
  if (button == GLUT_MIDDLE_BUTTON)
    {
      mygl.screen_voxel_pointed = NULL;
      return;
    }
  else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
    {
      mygl.screen_voxel_pointed = NULL;
      
      if (last_triangle.triangle_id != -1)
	{
	  display_id == 2;
	}
      last_triangle.boundary_id = -1;
      last_triangle.triangle_id = -1;
      return;
    }
  else if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
      for (n = 0; n < screen_voxels * screen_voxels; n++)
	{
	  mygl.screen_to_boundaries_map[ n ].vertex_z = 0.F;
	}
      for (boundary_index = 0; boundary_index < BOUNDARIES; boundary_index++)
	{
	  for (triangle_index = 0; triangle_index < mygl.boundary[ boundary_index ].triangles; triangle_index++)
	    {
	      triangle_p = &mygl.boundary[ boundary_index ].triangle[ triangle_index ];
	      
	      for (n = 0; n < 3; n++)
		{
		  myglTransformVertex (triangle_p->v[n].pos_x, triangle_p->v[n].pos_y, triangle_p->v[n].pos_z,
				       &x1, &y1, &z1);
		  
		  vertex_pixel_i = (int)(screen_voxels_screen_max_inv_x * (x1 + screen.max_x));
		  vertex_pixel_j = (int)(screen_voxels_screen_max_inv_y * (y1 + screen.max_y));
		  
		  if (vertex_pixel_i < 0 || vertex_pixel_i >= screen_voxels ||
		      vertex_pixel_j < 0 || vertex_pixel_j >= screen_voxels)
		    {
		      continue;
		    }
		  screen_voxel_p = &mygl.screen_to_boundaries_map[ vertex_pixel_i * screen_voxels + vertex_pixel_j ];
		  
		  if (z1 > screen_voxel_p->vertex_z)
		    {
		      screen_voxel_p->boundary_id = boundary_index;
		      screen_voxel_p->triangle_id = triangle_index;
		      screen_voxel_p->vertex_id   = n;
		      screen_voxel_p->vertex_z    = z1;
		    }
		}
	    }
	}
      mygl.screen_voxel_pointed = NULL;
      
      screen_voxel_p = &mygl.screen_to_boundaries_map[ mouse_pixel_i * screen_voxels + mouse_pixel_j ];
      
      if (screen_voxel_p->vertex_z > EPSILON)
	{
	  mygl.screen_voxel_pointed = screen_voxel_p;
	  
	  mygl.screen_voxel_coords.i = mouse_pixel_i;
	  mygl.screen_voxel_coords.j = mouse_pixel_j;
	  
	  passive_mouse_pixel_i = mouse_pixel_i;
	  passive_mouse_pixel_j = mouse_pixel_j;
	}
      else
	{
	  return;
	}
    }
  else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
      if (mygl.screen_voxel_pointed == NULL)
	{
	  if (display_id == 1)
	    {
	      display_id = 2;
	    }
	  return;
	}
      x1 = screen.max_x * (-1.F + (x0 << 1) / (float)viewport_pixels_x);
      y1 = screen.max_y * (-1.F + (y0 << 1) / (float)viewport_pixels_y);
      z1 = mygl.screen_voxel_pointed->vertex_z;
      
      myglAntiTransformVertex (x1, y1, z1, &x2, &y2, &z2);
      
      i = mygl.screen_voxel_pointed->vertex_id;
      j = mygl.screen_voxel_pointed->triangle_id;
      k = mygl.screen_voxel_pointed->boundary_id;
      
      triangle_p = &mygl.boundary[ k ].triangle[ j ];
      triangle_p->v[ i ].pos_x = x2;
      triangle_p->v[ i ].pos_y = y2;
      triangle_p->v[ i ].pos_z = z2;
      
      mygl.screen_voxel_coords.i = (int)(screen_voxels * (float)x0 / (float)viewport_pixels_x);
      mygl.screen_voxel_coords.j = (int)(screen_voxels * (float)y0 / (float)viewport_pixels_y);
      
      last_triangle.triangle_id = j;
      last_triangle.boundary_id = k;
      
      passive_mouse_pixel_i = mouse_pixel_i;
      passive_mouse_pixel_j = mouse_pixel_j;
    }
  else if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
    {
      if (last_triangle.triangle_id != -1)
	{
	  display_id == 2;
	}
      last_triangle.boundary_id = -1;
      last_triangle.triangle_id = -1;
      
      screen_voxel_p = &mygl.screen_to_boundaries_map[ mouse_pixel_i * screen_voxels + mouse_pixel_j ];
      
      if (screen_voxel_p->site_i < 0)
	{
	  return;
	}
      mygl.screen_voxel_pointed = screen_voxel_p;
      
      myglFromSiteToVoxelCoords (screen_voxel_p->site_i, screen_voxel_p->site_j, screen_voxel_p->site_k,
				 &mygl.selected_pixel_x, &mygl.selected_pixel_y, &mygl.selected_slice);
      
      myglFluidSitesIterativeSearching (mygl.selected_pixel_x,
					mygl.selected_pixel_y,
					mygl.selected_slice,
					mygl.selected_gray);
      
      myglSetBoundaryConfigurations ();
    }
  else
    {
      return;
    }
  display_id = 2;
}

void GLUTCALLBACK MotionFunction (int x0, int y0)
{
  float x1, y1, z1;
  float x2, y2, z2;
  
  int mouse_pixel_i, mouse_pixel_j;
  int i, j, k;
  
  Triangle *triangle_p;
  
  
  if (display_id == 1) return;
  
  y0 = viewport_pixels_y - y0;
  
  mouse_pixel_i = (int)(screen_voxels * (float)x0 / (float)viewport_pixels_x);
  mouse_pixel_j = (int)(screen_voxels * (float)y0 / (float)viewport_pixels_y);

  if (mouse_pixel_i < 0 || mouse_pixel_i >= screen_voxels ||
      mouse_pixel_j < 0 || mouse_pixel_j >= screen_voxels)
    {
      mygl.screen_voxel_pointed = NULL;
      
      //display_id = 0;
    }
  else if (mygl.screen_voxel_pointed != NULL &&
	   mygl.screen_voxel_pointed->vertex_z >= EPSILON)
    {
      x1 = screen.max_x * (-1.F + (x0 << 1) / (float)viewport_pixels_x);
      y1 = screen.max_y * (-1.F + (y0 << 1) / (float)viewport_pixels_y);
      z1 = mygl.screen_voxel_pointed->vertex_z;
      
      myglAntiTransformVertex (x1, y1, z1, &x2, &y2, &z2);
      
      i = mygl.screen_voxel_pointed->vertex_id;
      j = mygl.screen_voxel_pointed->triangle_id;
      k = mygl.screen_voxel_pointed->boundary_id;
      
      triangle_p = &mygl.boundary[ k ].triangle[ j ];
      triangle_p->v[ i ].pos_x = x2;
      triangle_p->v[ i ].pos_y = y2;
      triangle_p->v[ i ].pos_z = z2;
      
      mygl.screen_voxel_coords.i = (int)(screen_voxels * (float)x0 / (float)viewport_pixels_x);
      mygl.screen_voxel_coords.j = (int)(screen_voxels * (float)y0 / (float)viewport_pixels_y);
      
      last_triangle.triangle_id = j;
      last_triangle.boundary_id = k;
      
      display_id = 2;
    }
}


void GLUTCALLBACK PassiveMotionFunction (int x0, int y0)
{
  float x1, y1, z1;
  
  int mouse_pixel_i, mouse_pixel_j;
  int vertex_pixel_i, vertex_pixel_j;
  int n;
  int triangle_index;
  
  unsigned int boundary_index;
  
  Triangle *triangle_p;
  
  ScreenVoxel *screen_voxel_p;
  
  
  if (display_id == 1) return;
  
  //display_id = 0;
  
  passive_mouse_pixel_i = -1;
  passive_mouse_pixel_j = -1;
  
  y0 = viewport_pixels_y - y0;
  
  mouse_pixel_i = (int)(screen_voxels * (float)x0 / (float)viewport_pixels_x);
  mouse_pixel_j = (int)(screen_voxels * (float)y0 / (float)viewport_pixels_y);
  
  if (mouse_pixel_i < 0 || mouse_pixel_i >= screen_voxels ||
      mouse_pixel_j < 0 || mouse_pixel_j >= screen_voxels)
    {
      mygl.screen_voxel_pointed = NULL;
      return;
    }
  
  for (n = 0; n < screen_voxels * screen_voxels; n++)
    {
      mygl.screen_to_boundaries_map[ n ].vertex_z = 0.F;
    }
  for (boundary_index = 0; boundary_index < BOUNDARIES; boundary_index++)
    {
      for (triangle_index = 0; triangle_index < mygl.boundary[ boundary_index ].triangles; triangle_index++)
	{
	  triangle_p = &mygl.boundary[ boundary_index ].triangle[ triangle_index ];
	  
	  for (n = 0; n < 3; n++)
	    {
	      myglTransformVertex (triangle_p->v[n].pos_x, triangle_p->v[n].pos_y, triangle_p->v[n].pos_z,
				   &x1, &y1, &z1);
	      
	      vertex_pixel_i = (int)(screen_voxels_screen_max_inv_x * (x1 + screen.max_x));
	      vertex_pixel_j = (int)(screen_voxels_screen_max_inv_y * (y1 + screen.max_y));
	      
	      if (vertex_pixel_i < 0 || vertex_pixel_i >= screen_voxels ||
		  vertex_pixel_j < 0 || vertex_pixel_j >= screen_voxels)
		{
		  continue;
		}
	      screen_voxel_p = &mygl.screen_to_boundaries_map[ vertex_pixel_i * screen_voxels + vertex_pixel_j ];
	      
	      if (z1 > screen_voxel_p->vertex_z)
		{
		  screen_voxel_p->boundary_id = boundary_index;
		  screen_voxel_p->triangle_id = triangle_index;
		  screen_voxel_p->vertex_id   = n;
		  screen_voxel_p->vertex_z    = z1;
		}
	    }
	}
    }
  screen_voxel_p = &mygl.screen_to_boundaries_map[ mouse_pixel_i * screen_voxels + mouse_pixel_j ];
  
  if (screen_voxel_p->site_i >= 0)
    {
      passive_mouse_pixel_i = mouse_pixel_i;
      passive_mouse_pixel_j = mouse_pixel_j;
    }
  if (screen_voxel_p->vertex_z > EPSILON)
    {
      mygl.screen_voxel_coords.i = mouse_pixel_i;
      mygl.screen_voxel_coords.j = mouse_pixel_j;
      
      mygl.screen_voxel_pointed = screen_voxel_p;
      
      last_triangle.triangle_id = mygl.screen_voxel_pointed->triangle_id;
      last_triangle.boundary_id = mygl.screen_voxel_pointed->boundary_id;
    }
  else
    {
      return;
    }
  display_id = 2;
}


void GLUTCALLBACK Reshape (GLsizei w, GLsizei h)
{
  // the window is reshaped if necessary
  
  ortho_x *= (float)w / (float)viewport_pixels_x;
  ortho_y *= (float)h / (float)viewport_pixels_y;
  
  viewport_pixels_x = w;
  viewport_pixels_y = h;
  
  glViewport(0, 0, w, h);
  
  Projection ();
}


void myInit (void)
{
  mygl.screen_voxel_pointed = NULL;
  mygl.screen_voxel_coords.i = -1;
  
  last_triangle.triangle_id = -1;
  last_triangle.boundary_id = -1;
  
  ortho_x = 0.5F * mygl.system_size;
  ortho_y = 0.5F * mygl.system_size;
  viewpoint_radius = 10.F * mygl.system_size;
  
  Projection ();
}


void usage (char *progname)
{
  printf ("Usage: %s input path, output config, pars and checkpoint names\n", progname);
  printf ("slice and pixel size (mm), dataset resolution factor\n");
  printf ("    or\n");
  printf ("checkpoint file name, output config and pars file names\n");
}


int main (int argc, char *argv[])
{
  int is_checkpoint;
  
  
  if (argc == 8)
    {
      is_checkpoint = 0;
    }
  else if (argc == 4)
    {
      is_checkpoint = 1;
    }
  else
    {
      usage(argv[0]);
      exit(1);
    }
  
  if (!is_checkpoint)
    {
      input_path         = argv[1];
      output_config_name = argv[2];
      output_pars_name   = argv[3];
      checkpoint_name    = argv[4];
      
      slice_size = atof(argv[5]);
      pixel_size = atof(argv[6]);
      res_factor = atof(argv[7]);
      
      mygl.selected_slice = mygl.input_slices >> 1;
      mygl.selected_gray = 245;
      
      ReadConfig ();
      
      myglInitBoundaries ();
      
      display_id = 1;
    }
  else
    {
      checkpoint_name    = argv[1];
      output_config_name = argv[2];
      output_pars_name   = argv[3];
      
      ReadCheckpoint (checkpoint_name);
      
      display_id = 2;
    }
  
  
  glutInit (&argc, argv);
  
  myglOpenWindow (viewport_pixels_x, viewport_pixels_y);
  
  myInit ();
  
  
  glutReshapeFunc (Reshape);
  glutIdleFunc (Display);
  glutDisplayFunc (Display);
  glutKeyboardFunc (KeybordFunction);
  glutSpecialFunc (SpecialKeybordFunction);
  glutMouseFunc (MouseFunction);
  glutMotionFunc (MotionFunction);
  glutPassiveMotionFunc (PassiveMotionFunction);
  glutMainLoop ();

  return(0);
}
