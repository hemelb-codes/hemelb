/* This file contains the globals.
 * Should be removed ASAP
 */
#include "vis/colpixel.h"
#include "vis/rt.h"

#ifndef NO_STEER
#include <semaphore.h>
namespace heme
{
  namespace steering
  {
    pthread_mutex_t network_buffer_copy_lock;
    pthread_mutex_t LOCK;
    pthread_cond_t network_send_frame;
    sem_t nrl;
    sem_t connected_sem;
    sem_t steering_var_lock;
    
    bool is_frame_ready;
    bool connected;
    bool sending_frame;
    
    int send_array_length;
    
    int frame_size;
    
    bool updated_mouse_coords;
  }
}
#endif //NO_STEER


unsigned char *pixel_data = NULL;


double *f_old = NULL, *f_new = NULL;

int *f_id = NULL;






// TODO moving these requires making the collision / streaming functions non-static, which is potentially a big deal.
//double *inlet_density;
//double *outlet_density;
double* inlet_density, *outlet_density;


// 3 buffers needed for convergence-enabled simulations
double *f_to_send = NULL;
double *f_to_recv = NULL;

int *f_send_id = NULL;

int *f_recv_iv = NULL;

short int *f_data = NULL;

int is_bench;

// 3 variables needed for convergence-enabled simulations
double conv_error;
int cycle_tag, check_conv;
int is_inlet_normal_available;


int sites_x, sites_y, sites_z;
int blocks_x, blocks_y, blocks_z;
int blocks_yz, blocks;
int block_size, block_size2, block_size3, block_size_1;
int shift;
int sites_in_a_block;

double lbm_stress_type;
double lbm_stress_par;
double lbm_density_min, lbm_density_max;
double lbm_velocity_min, lbm_velocity_max;
double lbm_stress_min, lbm_stress_max;
double *lbm_average_inlet_velocity = NULL;
double *lbm_peak_inlet_velocity = NULL;
double *lbm_inlet_normal = NULL;
long int *lbm_inlet_count = NULL;

int lbm_terminate_simulation;

int net_machines;

namespace heme
{
  namespace vis
  {
    int col_pixels, col_pixels_max;
    int col_pixels_recv[2];
  
    int *col_pixel_id = NULL;
  
    ColPixel col_pixel_send[COLOURED_PIXELS_MAX];
    ColPixel col_pixel_recv[2][COLOURED_PIXELS_MAX];
  
    double pressure_min = 0.0, pressure_max = 0.0;
    double velocity_min = 0.0, velocity_max = 0.0;
    double stress_min = 0.0, stress_max = 0.0;
    double time = 0.0;

    int time_step = 0, cycle = 0;
    int period = 0, inlets = 0;
    int image_freq = 0;
    int pixels_max = 0;
    int streaklines = 1;

    float block_size_f;
    float block_size_inv;

    float physical_pressure_threshold_min;
    float physical_pressure_threshold_max;
    float physical_velocity_threshold_max;
    float physical_stress_threshold_max;
    float density_threshold_min, density_threshold_minmax_inv;
    float velocity_threshold_max_inv;
    float stress_threshold_max_inv;
    float brightness;
    float ctr_x, ctr_y, ctr_z;
    double mouse_pressure, mouse_stress;
    float streaklines_per_pulsatile_period, streakline_length;

    int mouse_x, mouse_y;
    int perform_rendering;
    int mode;


    Screen screen;

    Viewpoint viewpoint;

    Vis vis;
  }

}
