/* This file contains the globals.
 * Should be removed ASAP
 */
//#include "vis/colpixel.h"

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


double *f_old = 0, *f_new = 0;

int *f_id = 0;






// TODO moving these requires making the collision / streaming functions non-static, which is potentially a big deal.
//double *inlet_density;
//double *outlet_density;
double* inlet_density, *outlet_density;


// 3 buffers needed for convergence-enabled simulations
double *f_to_send = 0;
double *f_to_recv = 0;

int *f_send_id = 0;

int *f_recv_iv = 0;

short int *f_data = 0;

int is_bench;

// 3 variables needed for convergence-enabled simulations
double conv_error;
int cycle_tag, check_conv;
int is_inlet_normal_available;


int sites_x, sites_y, sites_z;
int blocks_x, blocks_y, blocks_z;
int blocks;
int block_size;
int shift;
int sites_in_a_block;

double lbm_stress_type;
double lbm_stress_par;
double lbm_density_min, lbm_density_max;
double lbm_velocity_min, lbm_velocity_max;
double lbm_stress_min, lbm_stress_max;
double *lbm_average_inlet_velocity = 0;
double *lbm_peak_inlet_velocity = 0;
double *lbm_inlet_normal = 0;
long int *lbm_inlet_count = 0;

int lbm_terminate_simulation;

int net_machines;

double lbm_phys_pressure_min = 0.0, lbm_phys_pressure_max = 0.0;
double lbm_phys_velocity_min = 0.0, lbm_phys_velocity_max = 0.0;
double lbm_phys_stress_min = 0.0, lbm_phys_stress_max = 0.0;

