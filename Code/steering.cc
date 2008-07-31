#include "config.h"
#include "network.h"
#include "steering.h"

float steer_par[STEERABLE_PARAMETERS+1] = {0.,0.,0.,    // scene center (dx,dy,dz)
					   45.0,45.0,   // longitude and latitude
					   1., 0.03,    // zoom and brightness
					   0.1, 0.1,    // velocity and stress ranges
					   -1., -1.,    // x-y position of the mouse of the client
					   0.,          // signal useful to terminate the simulation
					   0.};         // doRendering

void* hemeLB_steer (void* ptr)
{
  long int read_fd = (long int)ptr;
  
  printf("Kicking off steering thread with FD %i\n", (int)read_fd);
  
  while(1) {
    
    int num_chars = STEERABLE_PARAMETERS * sizeof(float) / sizeof(char);
    int bytes = sizeof(char) * num_chars;
    
    char* xdr_steering_data = (char*)malloc(bytes);
    
    XDR xdr_steering_stream;
    
    xdrmem_create(&xdr_steering_stream, xdr_steering_data, bytes, XDR_DECODE);
    
    int ret = recv_all(read_fd, xdr_steering_data, &num_chars);
    
    if (ret < 0) {
      // printf("Steering thread: broken network pipe...\n");
      break;
    }
    
    for (int i = 0; i < STEERABLE_PARAMETERS; i++)
      xdr_float(&xdr_steering_stream, &steer_par[i]);
    
    // printf("Got steering params ");
    // for (int i = 0; i < STEERABLE_PARAMETERS; i++) 
    //   printf("%0.4f ", steer_par[i]);
    // printf("\n"); 
    
    free(xdr_steering_data);
  }
  return 0;
}


void UpdateSteerableParameters (int *vis_perform_rendering, Vis *vis, LBM* lbm)
{
  steer_par[ STEERABLE_PARAMETERS ] = *vis_perform_rendering;
  
#ifndef NOMPI
    MPI_Bcast (steer_par, STEERABLE_PARAMETERS+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  float longitude, latitude;
  float zoom;
  float velocity_max, stress_max;
  float lattice_velocity_max, lattice_stress_max;
  
  
  vis_ctr_x     += steer_par[ 0 ];
  vis_ctr_y     += steer_par[ 1 ];
  vis_ctr_z     += steer_par[ 2 ];
  longitude      = steer_par[ 3 ];
  latitude       = steer_par[ 4 ];
  zoom           = steer_par[ 5 ];
  vis_brightness = steer_par[ 6 ];
  velocity_max   = steer_par[ 7 ];
  stress_max     = steer_par[ 8 ];
  
  vis_mouse_x              = (int)steer_par[  9 ];
  vis_mouse_y              = (int)steer_par[ 10 ];
  lbm_terminate_simulation = (int)steer_par[ 11 ];
  *vis_perform_rendering   = (int)steer_par[ 12 ];
  
  visConvertThresholds (velocity_max, stress_max,
			&lattice_velocity_max, &lattice_stress_max, lbm);
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
  		 PIXELS_X, PIXELS_Y,
  		 vis_ctr_x, vis_ctr_y, vis_ctr_z,
  		 5.F * vis->system_size,
  		 longitude, latitude,
  		 0.5F * (5.F * vis->system_size),
  		 zoom);
  
  vis_velocity_threshold_max_inv = 1.0/lattice_velocity_max;
  vis_stress_threshold_max_inv   = 1.0/lattice_stress_max;
}

