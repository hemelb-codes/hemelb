// In this file, the functions useful to initiate/end the LB simulation
// and perform the dynamics are reported

#include "config.h"

RT rt;

#include <string.h>

#ifdef RG

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <sys/wait.h>
#include <signal.h>


#define MYPORT 65250
#define CONNECTION_BACKLOG 10

#endif // RG

#ifdef RG
char host_name[255];

void *hemeLB_network (void *ptr)
{

	gethostname(host_name, 255);

#ifndef STEER
	FILE *f;
	f = fopen("env_details.asc","w");
	fprintf(f, "%s\n", host_name);
	fclose(f);

	printf("MPI 0 Hostname -> %s\n", host_name);
#endif

	signal(SIGPIPE, SIG_IGN); // Ignore a broken pipe
  
  int sock_fd;
  int new_fd;
  int yes = 1;
  
  int brokenPipe = 0;
  
  while(1) {
    
    pthread_mutex_lock ( &network_buffer_copy_lock );
    
    struct sockaddr_in my_address;
    struct sockaddr_in their_addr; // client address
    socklen_t sin_size;
    
    
    if ((sock_fd = socket (AF_INET, SOCK_STREAM, 0)) == -1) { perror("socket"); exit(1); }
    
    if (setsockopt (sock_fd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1) { perror("setsockopt"); exit(1); }
    
    my_address.sin_family = AF_INET;
    my_address.sin_port = htons (MYPORT);
    my_address.sin_addr.s_addr = INADDR_ANY;
    memset (my_address.sin_zero, '\0', sizeof my_address.sin_zero);
    
    if (bind (sock_fd, (struct sockaddr *)&my_address, sizeof my_address) == -1) { perror("bind"); exit(1); }
    
    if (listen(sock_fd, CONNECTION_BACKLOG) == -1) { perror("listen"); exit(1); }
    
    sin_size = sizeof their_addr;
    
    if ((new_fd = accept (sock_fd, (struct sockaddr *)&their_addr, &sin_size)) == -1) { perror("accept"); continue; }
    
    printf("server: got connection from %s\n", inet_ntoa (their_addr.sin_addr));
    
    close(sock_fd);
    
    int frameNumber = 0;
    
    brokenPipe = 0;
    
    pthread_mutex_unlock ( &network_buffer_copy_lock );
    
    while ( brokenPipe == 0 ) {
      
      pthread_mutex_lock ( &network_buffer_copy_lock );
      pthread_cond_wait (&network_send_frame, &network_buffer_copy_lock);
      
      int bytesSent = 0;
      
      XDR xdr_network_stream;
      
      u_int sizeToSend = 1024 * 1024;
      
      char* xdrBuffer = (char*) malloc( sizeToSend );
      
      xdrmem_create(&xdr_network_stream, xdrBuffer, sizeToSend, XDR_ENCODE);
      
      xdr_int(&xdr_network_stream, &frameNumber);
      xdr_int(&xdr_network_stream, &compressed_frame_size);
      
      for (int i = 0; i < compressed_frame_size; i++)
	{
	  xdr_u_char(&xdr_network_stream, &compressed_data[i]);
	}
      int currentPosition = xdr_getpos(&xdr_network_stream);
      
      int nElements = currentPosition / sizeof(char);
      
      printf("n Elements of char to send %i\n", nElements);
      
      for (int i = 0; i < nElements ; i++)
	{
	  int ret = send(new_fd, &xdrBuffer[i], sizeof(char), 0);
	  
	  if( ret < 0 ) {
	    
	    brokenPipe = 1;
	    break;
	    
	  } else {
	    
	    bytesSent += ret;
	  }
	}
      
      printf("bytes sent.... %i %i\n", bytesSent, nElements);
      
      xdr_destroy(&xdr_network_stream);
      
      free(xdrBuffer);
      
      printf("done sending array...");
      
      pthread_mutex_unlock ( &network_buffer_copy_lock );
      
      frameNumber++;
      
    } // while( brokenPipe == 0 )
    
    close(new_fd);
    
  } // while(1)
}

#endif // RG


inline void AbsorptionCoefficients (float flow_field_value, float t1, float t2, float cutoff, float *r, float *g, float *b)
{
  // the absorption factors are regulated here
  
  flow_field_value = fminf(1.F, flow_field_value);
  
  if (rt.is_isosurface)
    {
      *r += (1.F - flow_field_value) * (1.F - flow_field_value);
      *g += flow_field_value * (1.F - flow_field_value);
      *b += flow_field_value * flow_field_value;
    }
  else
    {
      // volume rendering option;
      // dt is the thickness of the trasversal segment in lattice unit
      
      float dt = t2 - t1;
      
      if (flow_field_value > cutoff)
	{
	  *r += dt * (1.F - flow_field_value) * (1.F - flow_field_value);
	  *g += dt * (flow_field_value) * (1.F - flow_field_value);
	  *b += dt * (flow_field_value * flow_field_value);
	}
      else
	{
	  dt *= flow_field_value;
	  *r += dt;
	  *g += dt;
	  *b += dt;
	}
    }
}


void usage (char *progname)
{
  printf ("Usage: %s /path/to/input/files\n", progname);
  printf (" files which must be present: config.dat input.asc, pars.asc, rt_pars.asc\n");
}


int main (int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output
  
  double simulation_time, total_time;
  
  int time_step, stability, is_converged;
  int write_checkpoint, check_convergence, perform_rt;
  int n;
  int checkpoint_count = 0;
  int convergence_count = 0;
  int ray_tracing_count = 0;
  int required_args = 2;

#ifdef RG
  pthread_t network_thread;
  pthread_attr_t pthread_attrib;
#endif // RG

#ifdef STEER
  int    reg_num_cmds;
  int    reg_cmds[REG_INITIAL_NUM_CMDS];
  char** steer_changed_param_labels;
  char** steer_recvd_cmd_params;

  SteerParams steer;

  steer_changed_param_labels = Alloc_string_array(REG_MAX_STRING_LENGTH,
						REG_MAX_NUM_STR_PARAMS);
  steer_recvd_cmd_params = Alloc_string_array(REG_MAX_STRING_LENGTH,
					    REG_MAX_NUM_STR_CMDS);

  int reg_finished;
#endif // STEER

  LBM lbm;
  
  Net net;
  
  net.err = MPI_Init (&argc, &argv);
  net.err = MPI_Comm_size (MPI_COMM_WORLD, &net.procs);
  net.err = MPI_Comm_rank (MPI_COMM_WORLD, &net.id);

  if (argc != required_args)
    {
      if (net.id == 0) usage(argv[0]);
      
      net.err = MPI_Abort (MPI_COMM_WORLD, 1);
      net.err = MPI_Finalize ();
    }
 
  char *input_file_path( argv[1] );

  char input_config_name[256];
  char input_parameters_name[256];
  char output_config_name[256];
  char checkpoint_config_name[256];
  char rt_parameters_name[256];
  char output_image_name[256];
  
  strcpy ( input_config_name , input_file_path );
  strcat ( input_config_name , "/config.dat" );

  strcpy ( input_parameters_name , input_file_path );
  strcat ( input_parameters_name , "/pars.asc" );

  strcpy ( output_config_name , input_file_path );
  strcat ( output_config_name , "/out.dat" );

  strcpy ( checkpoint_config_name , input_file_path );
  strcat ( checkpoint_config_name , "/check.dat" );

  strcpy ( rt_parameters_name , input_file_path );
  strcat ( rt_parameters_name , "/rt_pars.asc" );

  strcpy ( output_image_name , input_file_path );
  strcat ( output_image_name , "/image.dat" );

  if (net.id == 0)
    {
      printf ("\n");
      printf ("***********************************************************\n");
      printf("Opening parameters file: %s\n", input_parameters_name);
      printf("Opening config file: %s\n", input_config_name);
      printf("Opening rt parameters file: %s\n", rt_parameters_name);
    }
  
  lbm.inlet_density = NULL;
  lbm.outlet_density = NULL;
  
#ifdef STEER
  lbmReadParameters (&steer, input_parameters_name, &lbm, &net);
#else
  lbmReadParameters (input_parameters_name, &lbm, &net);
#endif

#ifdef STEER
  // create the derived datatype for the MPI_Bcast
  int count = 27;
  int blocklengths[27] = {1, 1, REG_MAX_NUM_STR_CMDS, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype types[27] = {MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_DOUBLE, MPI_DOUBLE, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL, MPI_INTEGER, MPI_UB};

  // calculate displacements
  MPI_Aint disps[27];
  disps[0] = 0;
  for(int i = 1; i < count; i++) {
    switch(types[i - 1]) {
    case MPI_INTEGER:
      disps[i] = disps[i - 1] + (sizeof(int) * blocklengths[i - 1]);
      break;
    case MPI_DOUBLE:
      disps[i] = disps[i - 1] + (sizeof(double) * blocklengths[i - 1]);
      break;
    case MPI_REAL:
      disps[i] = disps[i - 1] + (sizeof(float) * blocklengths[i - 1]);
      break;
    }
  }
  
  MPI_Datatype MPI_steer_type;
  MPI_Type_struct(count, blocklengths, disps, types, &MPI_steer_type);
  MPI_Type_commit(&MPI_steer_type);

  // initialize the steering library
  if(net.id == 0) {
    Steering_enable(REG_TRUE);

    reg_num_cmds = 2;
    reg_cmds[0] = REG_STR_STOP;
    reg_cmds[1] = REG_STR_PAUSE_INTERNAL;
    steer.status = Steering_initialize("HemeLB", reg_num_cmds, reg_cmds);
  }

  // broadcast/collect status
  net.err = MPI_Bcast(&steer.status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

  // if broken, quit
  if(steer.status == REG_FAILURE) {
    net.err = MPI_Finalize();
    return REG_FAILURE;
  }

#endif // STEER

#ifdef RG
  
  send_frame_count = 0;
  
  if(net.id == 0)
    {
      pthread_mutex_init (&network_buffer_copy_lock, NULL);
      pthread_cond_init (&network_send_frame, NULL);
      
      pthread_attr_init (&pthread_attrib);
      pthread_attr_setdetachstate (&pthread_attrib, PTHREAD_CREATE_JOINABLE);
      
      pthread_create (&network_thread, &pthread_attrib, hemeLB_network, NULL);
    }

#endif // RG

  lbmInit (input_config_name, checkpoint_config_name, &lbm, &net);
  
  netInit (&lbm, &net, &rt);
  
#ifdef STEER
  rtReadParameters (&steer, rt_parameters_name, &rt, &net);
#else
  rtReadParameters (rt_parameters_name, &rt, &net);
#endif

  
  if (!lbm.is_checkpoint)
    {
      lbmSetOptimizedInitialConditions (&lbm, &net);
    }
  
  rtInit (output_image_name, &rt);
  
  for (n = 0; n < 8; n++) net.timing[ n ] = 0.;
  
  stability = STABLE;
  checkpoint_count = 0;
  convergence_count = 0;
  ray_tracing_count = 0;

  int rotate_view = 1;
  
#ifdef STEER

  steer.rotate = rotate_view;

  // register params with RealityGrid here
  if(net.id == 0) {
    // read only and only if displaying
#ifdef RG
    steer.status = Register_param("Display Host", REG_FALSE,
				  (void*)(&host_name), REG_CHAR, "", "");
#endif

    // LBM params
    steer.status = Register_param("Tau", REG_TRUE,
				  (void*)(&steer.tau), REG_DBL, "0.5", "");
    steer.status = Register_param("Tolerance", REG_TRUE,
				  (void*)(&steer.tolerance), REG_DBL, "0.0", "0.1");
    steer.status = Register_param("Max time steps",REG_TRUE,
				  (void*)(&steer.max_time_steps), REG_INT, "1", "");
    steer.status = Register_param("Convergence frequency", REG_TRUE,
				  (void*)(&steer.conv_freq), REG_INT, "1", "");
    steer.status = Register_param("Checkpoint frequency", REG_TRUE,
				  (void*)(&steer.check_freq), REG_INT, "1", "");
    
    // RT params
    steer.status = Register_param("Auto Rotate", REG_TRUE,
				  (void*)(&steer.rotate), REG_INT, "0", "5");
    steer.status = Register_param("X pixel size", REG_TRUE,
				  (void*)(&steer.pixels_x), REG_INT, "0", "1024");
    steer.status = Register_param("Y pixel size", REG_TRUE,
				  (void*)(&steer.pixels_y), REG_INT, "0", "1024");
    steer.status = Register_param("Longitude", REG_TRUE,
				  (void *)(&steer.longitude), REG_FLOAT, "", "");
    steer.status = Register_param("Latitude", REG_TRUE,
				  (void *)(&steer.latitude), REG_FLOAT, "", "");
    steer.status = Register_param("Zoom", REG_TRUE,
				  (void *)(&steer.zoom), REG_FLOAT, "0.0", "");
    steer.status = Register_param("Image output frequency", REG_TRUE,
				  (void*)(&steer.image_freq), REG_INT, "1", "");
    steer.status = Register_param("Flow field type", REG_TRUE,
				  (void*)(&steer.flow_field_type), REG_INT, "0", "2");
    steer.status = Register_param("Is isosurface", REG_TRUE,
				  (void*)(&steer.is_isosurface), REG_INT, "0", "1");
    steer.status = Register_param("Absorption factor", REG_TRUE,
				  (void *)(&steer.abs_factor), REG_FLOAT, "0.0", "");
    steer.status = Register_param("Cutoff", REG_TRUE,
				  (void *)(&steer.cutoff), REG_FLOAT, "0.0", "1.0");
    steer.status = Register_param("Max density", REG_TRUE,
				  (void *)(&steer.max_density), REG_FLOAT, "0.0", "");
    steer.status = Register_param("Max velocity", REG_TRUE,
				  (void *)(&steer.max_velocity), REG_FLOAT, "0.0", "");
    steer.status = Register_param("Max stress", REG_TRUE,
				  (void *)(&steer.max_stress), REG_FLOAT, "0.0", "");
  }

  // broadcast/collect status
  net.err = MPI_Bcast(&steer.status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

  // if broken, quit
  if(steer.status == REG_FAILURE) {
    net.err = MPI_Finalize();
    return REG_FAILURE;
  }

  reg_finished = 0;

  if(net.id == 0) {
    printf("STEER: RealityGrid library initialized and parameters registered.\n");
    fflush(stdout);
  }

#endif // STEER

  for (time_step = 1; time_step <= lbm.time_steps_max; time_step++)
    {
      write_checkpoint = 0;
      check_convergence = 0;
      perform_rt = 0;
      
      if(net.id == 0)
      	{
      	  printf("time step: %i\n", time_step);
      	  fflush (stdout);
      	}
      
#ifdef STEER
      // call steering control
      if(net.id == 0) {
	steer.status = Steering_control(time_step,
					&steer.num_params_changed,
					steer_changed_param_labels,
					&steer.num_recvd_cmds,
					steer.recvd_cmds,
					steer_recvd_cmd_params);
      }
      
      // broadcast/collect everything
      net.err = MPI_Bcast(&steer, 1, MPI_steer_type, 0, MPI_COMM_WORLD);

      if(steer.status != REG_SUCCESS) {
	printf("STEER: I am %d and I detected that Steering_control failed.\n", net.id);
	fflush(stdout);
	continue;
      }
            
      // process commands received
      for(int i = 0; i < steer.num_recvd_cmds; i++) {
	switch(steer.recvd_cmds[i]) {
	case REG_STR_STOP:
	  printf("STEER: I am %d and I've been told to STOP.\n", net.id);
	  fflush(stdout);
	  reg_finished = 1;
	  break;
	}
      } // end of command processing
      
      // process changed params
      // not bothered what changed, just copy across...
      
      if(steer.num_params_changed > 0) {
	printf("STEER: I am %d and I was told that %d params changed.\n", net.id, steer.num_params_changed);
	fflush(stdout);
	
	// update lbm params
	lbm.tau = steer.tau;
	lbm.tolerance = steer.tolerance;
	lbm.time_steps_max = steer.max_time_steps;
	lbm.convergence_frequency = steer.conv_freq;
	lbm.checkpoint_frequency  = steer.check_freq;
	lbm.viscosity = ((2.0 * lbm.tau - 1.0) / 6.0);
	lbm.omega = -1.0 / lbm.tau;
	lbm.stress_par = (1.0 - 1.0 / (2.0 * lbm.tau)) / sqrt(2.0);

	// update rt params
	rtProjection (0.5F * rt.system_size, 0.5F * rt.system_size,
		      steer.pixels_x, steer.pixels_y,
		      steer.ctr_x, steer.ctr_y, steer.ctr_z,
		      5.F * rt.system_size,
		      steer.longitude, steer.latitude,
		      0.5F * (5.F * rt.system_size),
		      steer.zoom);

	rt.image_frequency = steer.image_freq;
	rt.flow_field_type = steer.flow_field_type;
	rt.is_isosurface = steer.is_isosurface;
	rt.absorption_factor = steer.abs_factor;
	rt.cutoff = steer.cutoff;
	rt.flow_field_value_max_inv[ DENSITY  ] = 1.F / steer.max_density;
	rt.flow_field_value_max_inv[ VELOCITY ] = 1.F / steer.max_velocity;
	rt.flow_field_value_max_inv[ STRESS   ] = 1.F / steer.max_stress;

	rotate_view = steer.rotate;
      }
      // end of param processing

#endif // STEER
      
      if (++checkpoint_count >= lbm.checkpoint_frequency)
	{
	  write_checkpoint = 1;
	  checkpoint_count = 0;
	}
      if (++convergence_count >= lbm.convergence_frequency)
	{
	  check_convergence = 1;
	  convergence_count = 0;
	}
      if (++ray_tracing_count >= rt.image_frequency)
	{
	  perform_rt = 1;
	  ray_tracing_count = 0;
	}

      // do we want to rotate?
      if(perform_rt == 1 && rotate_view > 0) {
	steer.longitude += rotate_view;
	rtProjection (0.5F * rt.system_size, 0.5F * rt.system_size,
		      steer.pixels_x, steer.pixels_y,
		      steer.ctr_x, steer.ctr_y, steer.ctr_z,
		      5.F * rt.system_size,
		      steer.longitude, steer.latitude,
		      0.5F * (5.F * rt.system_size),
		      steer.zoom);
      }

      // Between the rtRayTracingA/B calls, do not change any ray tracing
      // parameters.

	int flag = 0; // pthread_mutex_trylock( &network_buffer_copy_lock );

     // if (perform_rt)
      if (flag==0 && perform_rt==1)
	{
	  rtRayTracingA (AbsorptionCoefficients, &net, &rt);
	}

      stability = lbmCycle (write_checkpoint, check_convergence, perform_rt,
			    &is_converged, &lbm, &net);
      
      if (flag==0 && perform_rt==1)
	{
	  rtRayTracingB (AbsorptionCoefficients, &net, &rt);

//	      pthread_mutex_unlock (&network_buffer_copy_lock);

     // pthread_cond_signal (&network_send_frame);

	}
      
      if (stability == UNSTABLE || is_converged) break;

#ifdef STEER
      if(reg_finished == 1) break;
#endif // STEER
    }

  time_step = min(time_step, lbm.time_steps_max);
  
  simulation_time = 0.;
  
  for (n = 0; n < 5; n++)
    {
      simulation_time += net.timing[ n ];
    }
  
  if (net.id == 0)
    {
      if (stability == STABLE)
  	{
  	  if (!is_converged)
  	    {
  	      printf (" ATTENTION: SIMULATION NOT CONVERGED\n");
  	    }
  	}
      else
  	{
  	  printf (" ATTENTION: INSTABILITY CONDITION OCCURRED\n");
  	  printf (" AFTER %i time steps\n", time_step);
  	}
      printf (" fluid sites: %4i, MLSUPS: %.3f, time steps: %i \n procs: %i, machines: %i\n",
	      lbm.total_fluid_sites,
	      time_step * 1.E-6 * lbm.total_fluid_sites / simulation_time, time_step,
	      net.procs, net.machines);
      
      printf ("Opening output config file: %s\n", output_config_name);
      fflush (stdout);
    }
  
  lbmWriteConfig (stability, output_config_name, 0, &lbm, &net);
  
  net.err = MPI_Barrier (MPI_COMM_WORLD);
  
  if (net.id == 0)
    {
      printf ("timings results (average seconds per single time)\n");
      printf ("1)  rank\n");
      printf ("2)  collisions at interface-dependent lattice sites\n");
      printf ("3)  copying to the buffers to send/recv\n");
      printf ("4)  communicational time\n");
      printf ("5)  collisions + streaming from inner sites\n");
      printf ("6)  streaming from interface-dependent lattice sites\n");
      printf ("7)  all-reduce + convergence test + checkpoint\n");
      printf ("8)  intra-machine ray tracing \n");
      printf ("9)  inter-machine ray tracing \n");
      printf ("10) output (1 time only, seconds) \n");
      fflush (stdout);
    }
  total_time = 0.;
  
  for (n = 0; n < 8; n++) total_time += net.timing[ n ];
  
  for (n = 0; n < 8; n++)
    {
      if (net.timing[ n ] < 1.e-5 * total_time)
	{
	  net.timing[ n ] = 0.;
	}
      else
	{
	  if (n < 5)
	    {
	      net.timing[ n ] /= time_step;
	    }
	  else if (n == 5)
	    {
	      net.timing[ n ] /= net.convergence_count;
	    }
	  else if (n == 6 || n == 7)
	    {
	      net.timing[ n ] /= rt.ray_tracing_count;
	    }
	}
    }
  
  printf ("%i %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e 1) -> 10)\n",
	  net.id,
	  net.timing[0], net.timing[1], net.timing[2], net.timing[3],
	  net.timing[4], net.timing[5], net.timing[6], net.timing[7], net.timing[8]);
  fflush (stdout);
  
  rtEnd (&rt);
  lbmEnd (&lbm);
  netEnd (&net, &rt);

#ifdef RG
  pthread_join (network_thread, NULL);
#endif // RG

#ifdef STEER
  if(net.id == 0)
    {
      Steering_finalize();
      printf("STEER: Steering_finalize() called.\n");
      fflush(stdout);
    }
#endif // STEER

  net.err = MPI_Finalize ();
  
  return(0);
}

