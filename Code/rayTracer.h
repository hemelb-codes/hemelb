#ifndef __rayTracer_h_
#define __rayTracer_h_

#include "constants.h"
#include "visualisationControl.h"
#include "net.h"
// TODO this could probably be reduced to the net class and some visualisation class.
#include "rt.h"

// Some sort of coordinates.
struct BlockLocation
{
  short int i, j, k;
};

struct Cluster
{
  float minmax_x[2], minmax_y[2], minmax_z[2];
  
  float x[3];
  
  unsigned short int blocks_x, blocks_y, blocks_z;
  unsigned short int block_min[3];
};

struct AABB
{
  float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
};


void rtInit (Net *net);
void rtUpdateRayData (float *flow_field, float ray_t, float ray_segment, void (*ColourPalette) (float value, float col[]));
void rtRayTracing (void (*ColourPalette) (float value, float col[]));
void rtUpdateClusterVoxel (int i, float density, float velocity, float stress);
void rtEnd (void);

// TODO UGH get rid of them.
extern int cluster_blocks_vec[3];
extern int cluster_blocks_z, cluster_blocks_yz, cluster_blocks;
extern int clusters;
extern float block_size_f;
extern float block_size_inv;
extern float ray_dir[3];
extern float ray_inv[3];
extern float ray_vel_col[3];
extern float ray_stress_col[3];
extern float ray_length;
extern float ray_t_min;
extern float ray_density;
extern float ray_stress;




class rayTracer : public visualisationLayer
{

  public:




};

#endif //__rayTracer_h_