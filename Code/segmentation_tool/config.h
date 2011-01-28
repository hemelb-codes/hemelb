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

typedef long triangle_id;
typedef long site_index; // Typically in an array
typedef long voxel_id;
typedef long site_id;
typedef long block_id;
typedef unsigned long config_data;

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
    double sin_latitude, cos_latitude;
    double dist;
};

struct Site
{
    config_data cfg;
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

    triangle_id triangles;
};

struct Hit
{
    double pos[3];
    double t;

    triangle_id tri_id, previous_triangle_id;
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
    triangle_id *tri_ids;
    triangle_id triangles;
};

struct Mesh
{
    MeshTriangle *triangle;

    Voxel *voxel;

    double dim[3], half_dim[3];
    double voxel_size;

    int voxels[4];
    triangle_id triangles, triangles_max;

    double centre[3];
};

struct Coord
{
    site_id x[3];
    site_id iters;
};

struct ScreenVoxel
{
    double z[2];

    site_id site[3];
    triangle_id t_id;

    block_id b_id;
    voxel_id v_id;
};

struct Mouse
{
    triangle_id t_id;
    block_id b_id;
    voxel_id v_id;

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
    voxel_id input_voxels[3];
    voxel_id output_voxels[3];
    site_id sites[3];
    block_id blocks[3];
    site_id tot_sites;
    block_id tot_blocks;
    site_id stack_sites, stack_sites_max;
    int coords[COORD_BUFFERS];
    double seed_pos[3];
    site_index seed_site[3];
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
    double voxel_size;
    double stress_type;

    int res_factor;
    int screen_voxels;
    int mode;

    double segmentation_time, fps;
    Mesh mesh;

    Site *stack_site;
    Block *block;
    Coord *coord[COORD_BUFFERS];
    Boundary boundary[4];
    ScreenVoxel *screen_voxel;
    Mouse mouse;
    Menu menu;

    char *input_file;
    char *output_config;
    char *output_pars;
    char *checkpoint;
    char *output_coords;

};

extern unsigned long SOLID_TYPE;
extern unsigned long FLUID_TYPE;
extern unsigned long INLET_TYPE;
extern unsigned long OUTLET_TYPE;

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

extern unsigned long SITE_TYPE_MASK;
extern unsigned long BOUNDARY_CONFIG_MASK;
extern unsigned long BOUNDARY_DIR_MASK;
extern unsigned long BOUNDARY_ID_MASK;

extern int e[14 * 3];

extern int inv_dir[14];

extern Screen screen;

extern Viewpoint viewpoint;

extern Vis vis;

#endif // CONFIG
