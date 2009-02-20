#ifndef CONFIG
#define CONFIG


#include <fstream>
#include <vector>

using namespace std;

#ifndef GLUTCALLBACK
#define GLUTCALLBACK
#endif


#define DEG_TO_RAD    0.01745329

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
  float col[3];
  float ctr[3];
  float dim[2];
  float zoom;
  
  int pixels[2];
};


struct Viewpoint
{
  float pos[3];
  
  float sin_longitude;
  float cos_longitude;
  
  float sin_latitude;
  float cos_latitude;
  
  float dist;
};


struct Site
{
  unsigned int cfg;
  unsigned int label;
};


struct Block
{
  Site *site;
};


struct Vertex
{
  float pos[3];
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
  
  float pos[3];
  float nor[3];
  
  float pressure_avg, pressure_amp, pressure_phs;
  
  int normal_sign;
};


struct Boundary
{
  Triangle *triangle;
  
  int triangles;
};


struct Coord
{
  short int x[3];
  short int iters;
};


struct ScreenVoxel
{
  float z[2];
  
  short int t_id;
  short int site[3];
  
  char b_id;
  char v_id;
};


struct Mouse
{
  short int t_id;
  char b_id;
  char v_id;
  
  short int x[2], dy;
  short int state;
};


struct Menu
{
  int id;
  int option;
};


struct Vis
{
  float scale[3], inv_scale[3];
  
  int input_voxels[3];
  int output_voxels[3];
  int pixel_depth;
  
  int sites[3];
  int blocks[3];
  int tot_sites, tot_blocks;
  int stack_sites, stack_sites_max;
  int coords[COORD_BUFFERS];
  
  short int selected_voxel[3];
  
  float selected_grey;
  float grey_min, grey_max;
  
  int viewport_pixels[2];
  
  float background[3];
  
  float ortho[2];
  float longitude, latitude;
  float viewpoint_radius;
  float viewport_radius;
  float zoom;
  float scene_center[3];
  
  float dim[3], half_dim[3];
  float system_size;
  
  float slice_size, pixel_size;
  
  int res_factor;
  int screen_voxels;
  int mode;
  
  float segmentation_time, fps;  
  
  Site *stack_site;
  
  Block *block;
  
  // unsigned short int *voxel;
  signed short int *voxel;
  
  Coord *coord[COORD_BUFFERS];
  
  Boundary boundary[4];
  
  ScreenVoxel *screen_voxel;
  
  Mouse mouse;
  
  Menu menu;
  
  char *input_path;
  char *output_config;
  char *output_pars;
  char *checkpoint;
  
  vector<string> file_list;

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
