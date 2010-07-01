#include "config.h"

// parameters related to the lattice directions

int e_x[] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
int e_y[] = { 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1};
int e_z[] = { 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1};
int inv_dir[] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};

#ifndef NOMPI
MPI_Datatype MPI_col_pixel_type;
#endif


#ifndef NO_STEER
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
#endif


unsigned char *pixel_data = NULL;


double *f_old = NULL, *f_new = NULL;

int *f_id = NULL;



// 3 buffers needed for convergence-enabled simulations
double *f_to_send = NULL;
double *f_to_recv = NULL;

int *f_send_id = NULL;

int *f_recv_iv = NULL;

short int *f_data = NULL;

double *inlet_density = NULL;
double *inlet_density_avg = NULL, *inlet_density_amp = NULL, *inlet_density_phs = NULL;
double *outlet_density = NULL;
double *outlet_density_avg = NULL, *outlet_density_amp = NULL, *outlet_density_phs = NULL;


int col_pixels, col_pixels_max;
int col_pixels_recv[2];

int *col_pixel_id = NULL;

ColPixel col_pixel_send[COLOURED_PIXELS_MAX];
ColPixel col_pixel_recv[2][COLOURED_PIXELS_MAX];

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

double vis_pressure_min = 0.0, vis_pressure_max = 0.0;
double vis_velocity_min = 0.0, vis_velocity_max = 0.0;
double vis_stress_min = 0.0, vis_stress_max = 0.0;
double vis_time = 0.0;

int vis_time_step = 0, vis_cycle = 0;
int vis_period = 0, vis_inlets = 0;
int vis_image_freq = 0;
int vis_pixels_max = 0;
int vis_streaklines = 1;


float block_size_f;
float block_size_inv;
float vis_physical_pressure_threshold_min;
float vis_physical_pressure_threshold_max;
float vis_physical_velocity_threshold_max;
float vis_physical_stress_threshold_max;
float vis_density_threshold_min, vis_density_threshold_minmax_inv;
float vis_velocity_threshold_max_inv;
float vis_stress_threshold_max_inv;
float vis_brightness;
float vis_ctr_x, vis_ctr_y, vis_ctr_z;
double vis_mouse_pressure, vis_mouse_stress;
double vis_glyph_length;
float vis_streaklines_per_pulsatile_period, vis_streakline_length;

int vis_mouse_x, vis_mouse_y;
int vis_perform_rendering;
int vis_mode;

int cluster_blocks_vec[3];
int cluster_blocks_z, cluster_blocks_yz, cluster_blocks;


float ray_dir[3];
float ray_inv[3];
float ray_vel_col[3];
float ray_stress_col[3];
float ray_length;
float ray_t_min;
float ray_density;
float ray_stress;

int clusters;

Screen screen;

Viewpoint viewpoint;

Vis vis;



