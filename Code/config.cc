#include "config.h"


float EPSILON = 1.0e-30;

int STABLE   = 1;
int UNSTABLE = 0;

// the constants needed to define the configuration of the lattice
// sites follow

unsigned int SOLID_TYPE  = 0U;
unsigned int FLUID_TYPE  = 1U;
unsigned int INLET_TYPE  = 2U;
unsigned int OUTLET_TYPE = 3U;
unsigned int NULL_TYPE   = 4U;

unsigned int BOUNDARIES              = 4U;
unsigned int INLET_BOUNDARY          = 0U;
unsigned int OUTLET_BOUNDARY         = 1U;
unsigned int WALL_BOUNDARY           = 2U;
unsigned int CHARACTERISTIC_BOUNDARY = 3U;

unsigned int SITE_TYPE_BITS       = 2U;
unsigned int BOUNDARY_CONFIG_BITS = 14U;
unsigned int BOUNDARY_DIR_BITS    = 4U;
unsigned int BOUNDARY_ID_BITS     = 10U;

unsigned int BOUNDARY_CONFIG_SHIFT = 2U;   // SITE_TYPE_BITS;
unsigned int BOUNDARY_DIR_SHIFT    = 16U;  // BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
unsigned int BOUNDARY_ID_SHIFT     = 20U;  // BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

unsigned int SITE_TYPE_MASK       = ((1U <<  2U) - 1U);         // ((1U << SITE_TYPE_BITS) - 1U);
unsigned int BOUNDARY_CONFIG_MASK = ((1U << 14U) - 1U) << 2U;   // ((1U << BOUNDARY_CONFIG_BITS) - 1U) << BOUNDARY_CONFIG_SHIFT;
unsigned int BOUNDARY_DIR_MASK    = ((1U <<  4U) - 1U) << 16U;  //((1U << BOUNDARY_DIR_BITS) - 1U)    << BOUNDARY_DIR_SHIFT;
unsigned int BOUNDARY_ID_MASK     = ((1U << 10U) - 1U) << 20U;  // ((1U << BOUNDARY_ID_BITS) - 1U)     << BOUNDARY_ID_SHIFT
unsigned int CHARACTERISTIC_MASK  = 1U << 31U;

// square of the speed of sound

double Cs2 = 3.0 / 8.0;


// parameters related to the lattice directions

int e_x[] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
int e_y[] = { 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1};
int e_z[] = { 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1};
int inv_dir[] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};


#ifdef RG

pthread_mutex_t network_buffer_copy_lock;
pthread_cond_t network_send_frame;

int send_array_length;
int send_frame_count;

int compressed_frame_size;

unsigned char *pixel_data = NULL;
unsigned char *compressed_data = NULL;

#endif // RG


double *f_old = NULL, *f_new = NULL;

int *f_id = NULL;

Velocity *vel = NULL;

float *flow_field = NULL;

double *d = NULL;

double **nd_p = NULL;


short int f_data[4*SHARED_DISTRIBUTIONS_MAX];

double f_to_send[SHARED_DISTRIBUTIONS_MAX];
double f_to_recv[SHARED_DISTRIBUTIONS_MAX];

int f_send_id[SHARED_DISTRIBUTIONS_MAX];
int f_recv_iv[SHARED_DISTRIBUTIONS_MAX];

float pixel_color_to_send[4*1024*1024];
float pixel_color_to_recv[4*1024*1024];


Screen screen;

Viewpoint viewpoint;

Ray ray;


// some simple functions

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


double myClock ()
{
  return (double)clock () / (double)CLOCKS_PER_SEC;
}
