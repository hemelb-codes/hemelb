#include "config.h"

double PI = 2. * atan2(1., 0.);

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

#ifndef NOMPI
MPI_Datatype MPI_col_pixel_type;
#endif


#ifdef RG

pthread_mutex_t network_buffer_copy_lock;
pthread_cond_t network_send_frame;

int send_array_length;

int frame_size;

unsigned char *pixel_data = NULL;

#endif // RG


double *f_old = NULL, *f_new = NULL;

int *f_id = NULL;

#ifndef TD
double *vel = NULL;
#endif

float *flow_field = NULL;

#ifndef BENCH
double *d = NULL;

double **nd_p = NULL;
#endif

Cluster *cluster = NULL;


short int f_data[4*SHARED_DISTRIBUTIONS_MAX];

double f_to_send[SHARED_DISTRIBUTIONS_MAX];
double f_to_recv[SHARED_DISTRIBUTIONS_MAX];

int f_send_id[SHARED_DISTRIBUTIONS_MAX];
int f_recv_iv[SHARED_DISTRIBUTIONS_MAX];

float streamline_to_send[NEIGHBOUR_PROCS_MAX][VIS_VEC_SIZE*STREAMLINES_MAX];
float streamline_to_recv[NEIGHBOUR_PROCS_MAX][1+VIS_VEC_SIZE*STREAMLINES_MAX];

int streamlines_to_send[NEIGHBOUR_PROCS_MAX];
int streamlines_to_recv[NEIGHBOUR_PROCS_MAX];


int sites_x, sites_y, sites_z;
int blocks_x, blocks_y, blocks_z;
int blocks_yz, blocks;
int block_size, block_size2, block_size3, block_size_1;
int shift;
int sites_in_a_block;

float block_size_inv;


float ray_dir[VIS_VEC_SIZE];
float ray_inv[VIS_VEC_SIZE];
float ray_col[VIS_VEC_SIZE];


short int clusters;


Screen screen;

Viewpoint viewpoint;

Vis vis;


// some simple math functions

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
//#ifndef NOMPI
  struct timeval time_data;
  
  gettimeofday (&time_data, NULL);
  
  return (double)time_data.tv_sec + (double)time_data.tv_usec / 1.e6;
//#else
//  return MPI_Wtime();
//#endif
  
  //double time;
  //
  //int rc;
  //
  //struct timespec ts;
  //
  //
  //rc = clock_gettime (CLOCK_REALTIME, &ts);
  //
  //if (rc == 0)
  //  {
  //    time = (double)ts.tv_sec + (double)ts.tv_nsec / 1000000000.;
  //  }
  //else
  //  {
  //    fprintf(stderr, "ERROR: clock_gettime() failed\n");
  //    time = 0.;
  //  }
  //return time;
}


