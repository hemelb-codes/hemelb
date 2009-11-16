#ifndef CONFIG
#define CONFIG

#ifdef HEMELB_OSX
	#include <GLUT/glut.h>
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <GL/glut.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>
#include "math.h"

using namespace std;

#ifndef GLUTCALLBACK
#define GLUTCALLBACK
#endif


#define EPSILON       1.0e-30
#define DEG_TO_RAD    0.01745329

#define VON_MISES_STRESS   +1.0
#define SHEAR_STRESS       -1.0

#define BLOCK_SIZE        8
#define SITES_PER_BLOCK   (BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE)
#define SHIFT             3


#define SUCCESS   1

#define COORD_BUFFERS         3
#define COORD_BUFFER_SIZE_A   1000000
#define COORD_BUFFER_SIZE_B   1000000
#define COORD_BUFFER_SIZE_C   1000000


#define VoxelId(voxel,voxels) (voxel[2]*voxels[1]+voxel[1])*voxels[0]+voxel[0]
#define BlockId(b_id,blocks) (b_id[0]*blocks[1]+b_id[1])*blocks[2]+b_id[2]
#define SiteId(site,b_id) ((((site[0]-(b_id[0]<<SHIFT))<<SHIFT)+site[1]-(b_id[1]<<SHIFT))<<SHIFT)+site[2]-(b_id[2]<<SHIFT)


struct Screen
{
  double col[3];
  double ctr[3];
  double dim[2];
  double zoom;
  
  int pixels[2];
};


struct Viewpoint
{
  double pos[3];
  double sin_longitude, cos_longitude;
  double sin_latitude,  cos_latitude;
  double dist;
};


struct Site
{
  unsigned int cfg;
  unsigned int label;
#ifdef MESH
  int triangle_id;
#endif
};


struct Block
{
  Site *site;
};


struct Vertex
{
  double pos[3];
};


struct Disc
{
  double sin_longitude, cos_longitude;
  double sin_latitude, cos_latitude;
  double r2;
};


struct BoundaryTriangle
{
  Vertex v[3];
  
  Disc d;
  
  double pos[3];
  double nor[3];
  
  double pressure_avg, pressure_amp, pressure_phs;
  
  int normal_sign;
};


struct Boundary
{
  BoundaryTriangle *triangle;
  
  int triangles;
};


#ifdef MESH
struct Hit
{
  double pos[3];
  double t;
  
  int triangle_id, previous_triangle_id;
};


struct Ray
{
  double org[3], dir[3];
  
  double t_max, t_near, t_far;
};


struct MeshTriangle
{
  Vertex v[3];
  
  double nor[3];
};


struct Voxel
{
  int *triangle_id;
  int triangles;
};


struct Mesh
{
  MeshTriangle *triangle;
  
  Voxel *voxel;
  
  double dim[3], half_dim[3];
  double voxel_size;
  
  int voxels[4];
  int triangles, triangles_max;
};
#endif // MESH


struct Coord
{
  short int x[3];
  short int iters;
};


struct ScreenVoxel
{
  double z[2];
  
  short int site[3];
  short int t_id;
  
  char b_id;
  char v_id;
};


struct Mouse
{
  short int t_id;
  char b_id;
  char v_id;
  
  short int x[2];
  short int dy;
  short int state;
};


struct Menu
{
  int id;
  int option;
};


struct Vis
{
  int input_voxels[3];
  int output_voxels[3];
#ifndef MESH
  int pixel_depth;
#endif
  int sites[3];
  int blocks[3];
  int tot_sites, tot_blocks;
  int stack_sites, stack_sites_max;
  int coords[COORD_BUFFERS];
#ifndef MESH
  short int selected_voxel[3];
  
  double selected_grey;
  double grey_min, grey_max;
#else
  double seed_pos[3];
  
  short int seed_site[3];
#endif
  int viewport_pixels[2];
  
  double background[3];
  double ortho[2];
  double longitude, latitude;
  double viewpoint_radius;
  double viewport_radius;
  double zoom;
  double scene_center[3];
  
  double dim[3], half_dim[3];
  double system_size;
#ifndef MESH
  double slice_size, pixel_size;
#else
  double voxel_size;
#endif
  double stress_type;
  
  int res_factor;
  int screen_voxels;
  int mode;
  
  double segmentation_time, fps;  
#ifndef MESH
  short int *voxel;
#else
  Mesh mesh;
#endif
  Site *stack_site;
  
  Block *block;
  
  Coord *coord[COORD_BUFFERS];
  
  Boundary boundary[4];
  
  ScreenVoxel *screen_voxel;
  
  Mouse mouse;
  
  Menu menu;
  
#ifndef MESH
  char *input_path;
#else
  char *input_file;
#endif
  char *output_config;
  char *output_pars;
  char *checkpoint;
#ifndef MESH
  vector<string> file_list;
#endif

};


extern unsigned int SOLID_TYPE;
extern unsigned int FLUID_TYPE;
extern unsigned int INLET_TYPE;
extern unsigned int OUTLET_TYPE;

extern unsigned int BOUNDARIES;
extern unsigned int INLET_BOUNDARY;
extern unsigned int OUTLET_BOUNDARY;
extern unsigned int WALL_BOUNDARY;

extern unsigned int SITE_TYPE_BITS;
extern unsigned int BOUNDARY_CONFIG_BITS;
extern unsigned int BOUNDARY_DIR_BITS;
extern unsigned int BOUNDARY_ID_BITS;

extern unsigned int BOUNDARY_CONFIG_SHIFT;
extern unsigned int BOUNDARY_DIR_SHIFT;
extern unsigned int BOUNDARY_ID_SHIFT;

extern unsigned int SITE_TYPE_MASK;
extern unsigned int BOUNDARY_CONFIG_MASK;
extern unsigned int BOUNDARY_DIR_MASK;
extern unsigned int BOUNDARY_ID_MASK;


extern int e[14*3];

extern int inv_dir[14];


extern Screen screen;

extern Viewpoint viewpoint;

extern Vis vis;

#endif // CONFIG
