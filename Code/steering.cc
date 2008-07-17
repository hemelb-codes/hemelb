#include "config.h"
#include "network.h"
#include "steering.h"

float steer_par[STEERABLE_PARAMETERS] = {0.0, 0.0, 0.0, 45.0, 45.0, 1.0, 0.03, 0.001, 0.01};

void* hemeLB_steer (void* ptr) {

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

/*		printf("Got steering params ");
		for (int i = 0; i < STEERABLE_PARAMETERS; i++) 
			printf("%0.4f ", steer_par[i]);
		printf("\n"); */

		free(xdr_steering_data);
	}

	return 0;

}

void UpdateSteerableParameters(Vis *vis, LBM* lbm) {

#ifndef NOMPI
  MPI_Bcast (steer_par, STEERABLE_PARAMETERS, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  float ctr_x, ctr_y, ctr_z;
  float longitude, latitude;
  float zoom;
  float velocity_max, stress_max;
  float lattice_velocity_max, lattice_stress_max;
  
  ctr_x          = steer_par[ 0 ];
  ctr_y          = steer_par[ 1 ];
  ctr_z          = steer_par[ 2 ];
  longitude      = steer_par[ 3 ];
  latitude       = steer_par[ 4 ];
  zoom           = steer_par[ 5 ];
  vis_brightness = steer_par[ 6 ];
  velocity_max   = steer_par[ 7 ];
  stress_max     = steer_par[ 8 ];

  visConvertThresholds (velocity_max, stress_max,
			&lattice_velocity_max, &lattice_stress_max, lbm);
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
  		 PIXELS_X, PIXELS_Y,
  		 ctr_x, ctr_y, ctr_z,
  		 5.F * vis->system_size,
  		 longitude, latitude,
  		 0.5F * (5.F * vis->system_size),
  		 zoom);
  
  vis_velocity_threshold_max_inv = 1.0/velocity_max;
  vis_stress_threshold_max_inv   = 1.0/stress_max;

}

