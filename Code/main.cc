// In this file, the functions useful to initiate/end the LB simulation
// and perform the dynamics are reported

#include "config.h"


#ifdef RG
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#define MYPORT 65250
#define CONNECTION_BACKLOG 10
#endif // RG


FILE *timings_ptr;


#ifdef RG

char host_name[255];

u_int sizeToSend = (2 + SCREEN_SIZE_MAX) * sizeof(unsigned int);

char *xdrSendBuffer;


int send_all(int sockid, char *buf, int *length ) {
	
  int sent_bytes = 0;
  int bytes_left_to_send = *length;
  int n;
	
  while( sent_bytes < *length ) {
    n = send(sockid, buf+sent_bytes, bytes_left_to_send, 0);
    if (n == -1)
      break;
    sent_bytes += n;
    bytes_left_to_send -= n;
  }
	
  *length = sent_bytes;
	
  return n==-1?-1:0;

}

void *hemeLB_network (void *ptr)
{
  gethostname (host_name, 255);

#ifndef STEER
  FILE *f = fopen ("env_details.asc","w");
  
  fprintf (f, "%s\n", host_name);
  fclose (f);
  
  fprintf (timings_ptr, "MPI 0 Hostname -> %s\n\n", host_name);
#endif

  signal(SIGPIPE, SIG_IGN); // Ignore a broken pipe
  
  int sock_fd;
  int new_fd;
  int yes = 1;
  
  int is_broken_pipe = 0;
  int frame_number = 0;
  
  while (1)
    {
      pthread_mutex_lock ( &network_buffer_copy_lock );
      
      struct sockaddr_in my_address;
      struct sockaddr_in their_addr; // client address
      socklen_t sin_size;
      
      
      if ((sock_fd = socket (AF_INET, SOCK_STREAM, 0)) == -1)
	{
	  perror("socket");
	  exit (1);
	}
      
      if (setsockopt (sock_fd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1)
	{
	  perror("setsockopt");
	  exit (1);
	}
      
      my_address.sin_family = AF_INET;
      my_address.sin_port = htons (MYPORT);
      my_address.sin_addr.s_addr = INADDR_ANY;
      memset (my_address.sin_zero, '\0', sizeof my_address.sin_zero);
      
      if (bind (sock_fd, (struct sockaddr *)&my_address, sizeof my_address) == -1)
	{
	  perror ("bind");
	  exit (1);
	}
      
      if (listen (sock_fd, CONNECTION_BACKLOG) == -1)
	{
	  perror ("listen");
	  exit (1);
	}
      
      sin_size = sizeof their_addr;
      
      if ((new_fd = accept (sock_fd, (struct sockaddr *)&their_addr, &sin_size)) == -1)
	{
	  perror("accept");
	  continue;
	}
      
      fprintf (timings_ptr, "server: got connection from %s\n", inet_ntoa (their_addr.sin_addr));
      
      close(sock_fd);
      
      is_broken_pipe = 0;
      
      pthread_mutex_unlock ( &network_buffer_copy_lock );
      
      while (is_broken_pipe == 0)
	{
	  pthread_mutex_lock ( &network_buffer_copy_lock );
	  pthread_cond_wait (&network_send_frame, &network_buffer_copy_lock);
	  
	  int bytesSent = 0;
	  
	  XDR xdr_network_stream;
	  
	  unsigned int four_compressed_data, data;
	  
	  int m;
	  
	  
	  xdrmem_create (&xdr_network_stream, xdrSendBuffer, sizeToSend, XDR_ENCODE);
	  
	  xdr_int (&xdr_network_stream, &frame_number);
	  xdr_int (&xdr_network_stream, &compressed_frame_size);
	  
	  for (int i = 0; i < compressed_frame_size; i++)
	    {
	      xdr_u_char (&xdr_network_stream, &compressed_data[i]);
	    }
	  //m = (compressed_frame_size >> 2) << 2;
	  //
	  //for (int i = 0; i < m;)
	  //  {
	  //    four_compressed_data = compressed_data[i++];
	  //    four_compressed_data |= (data = (unsigned int)compressed_data[i++]) << 8U;
	  //    four_compressed_data |= (data = (unsigned int)compressed_data[i++]) << 16U;
	  //    four_compressed_data |= (data = (unsigned int)compressed_data[i++]) << 24U;
	  //    
	  //    xdr_u_int (&xdr_network_stream, &four_compressed_data);
	  //  }
	  //for (int i = m; i < compressed_frame_size;)
	  //  {
	  //    four_compressed_data = compressed_data[i++];
	  //    
	  //    if (i++ < compressed_frame_size)
	  //	{
	  //	  four_compressed_data |= (data = (unsigned int)compressed_data[i]) << 8U;
	  //	}
	  //    if (i++ < compressed_frame_size)
	  //	{
	  //	  four_compressed_data |= (data = (unsigned int)compressed_data[i]) << 16U;
	  //	}
	  //    if (i++ < compressed_frame_size)
	  //	{
	  //	  four_compressed_data |= (data = (unsigned int)compressed_data[i]) << 24U;
	  //	}
	  //    xdr_u_int (&xdr_network_stream, &four_compressed_data);
	  //  }

	  int currentPosition = xdr_getpos (&xdr_network_stream);
	  
	  int nElements = currentPosition / sizeof(char);
	  
	  fprintf (timings_ptr, "n Elements of char to send %i (%iB)\n", nElements, currentPosition);

          int ret;

          ret = send_all(new_fd, xdrSendBuffer, &currentPosition);

          if (ret < 0) {
            is_broken_pipe = 1;
            break;
          } else {
            bytesSent += currentPosition;
          }

	  fprintf (timings_ptr, "bytes sent, elements: %i %i\n", bytesSent, nElements);
	  
	  xdr_destroy (&xdr_network_stream);
	  
	  pthread_mutex_unlock ( &network_buffer_copy_lock );
	  
	  frame_number++;
	  
	} // while (is_broken_pipe == 0)
      
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
#ifndef BENCH
  fprintf (timings_ptr, "Usage: %s path of the input files\n", progname);
  fprintf (timings_ptr, "the following files must be present: config.dat, pars.asc\n");
  fprintf (timings_ptr, "check.dat (if it is a checkpoint file and must be used), rt_pars.asc\n");
#else
  fprintf (timings_ptr, "Usage: %s path of the input files and minutes for benchmarking\n", progname);
  fprintf (timings_ptr, "the following files must be present: config.dat, pars.asc\n");
  fprintf (timings_ptr, "check.dat (if it is a checkpoint file and must be used), rt_pars.asc\n");
#endif
}


int main (int argc, char *argv[])
{
  double total_time = myClock ();
  
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output
  
#ifndef BENCH
  double simulation_time;
#else
  double minutes;
  double fluid_solver_time, fluid_solver_and_vr_time, fluid_solver_and_is_time;
#endif // BENCH
  
  int time_step, stability, is_converged;
  int write_checkpoint, check_convergence, perform_rt;
  int n;
  int checkpoint_count = 0;
  int convergence_count = 0;
  int ray_tracing_count = 0;
  int is_bench_section_finished;
  int is_thread_locked;
  
#ifdef BENCH
  int fluid_solver_time_steps;
  int fluid_solver_and_vr_time_steps;
  int fluid_solver_and_is_time_steps;
#endif
  
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

  steer_changed_param_labels = Alloc_string_array (REG_MAX_STRING_LENGTH,
						   REG_MAX_NUM_STR_PARAMS);
  steer_recvd_cmd_params = Alloc_string_array (REG_MAX_STRING_LENGTH,
					       REG_MAX_NUM_STR_CMDS);
  
  int reg_finished;
#endif // STEER

  LBM lbm;
  
  Net net;
  
  
  net.err = MPI_Init (&argc, &argv);
  net.err = MPI_Comm_size (MPI_COMM_WORLD, &net.procs);
  net.err = MPI_Comm_rank (MPI_COMM_WORLD, &net.id);
  
  if (argc != 2 && argc != 3)
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
  char timings_name[256];
  
  
#ifdef BENCH
  minutes = atof( argv[2] );
#endif
  
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
  
  strcpy ( timings_name , input_file_path );
  strcat ( timings_name , "/timings.asc" );
  
  
  if (net.id == 0)
    {
      timings_ptr = fopen (timings_name, "w");
      
      fprintf (timings_ptr, "***********************************************************\n");
      fprintf (timings_ptr, "Opening parameters file:\n %s\n", input_parameters_name);
      fprintf (timings_ptr, "Opening config file:\n %s\n", input_config_name);
      fprintf (timings_ptr, "Opening rt parameters file:\n %s\n\n", rt_parameters_name);
    }
  
  lbm.inlet_density = NULL;
  lbm.outlet_density = NULL;
  
#ifdef STEER
  lbmReadParameters (input_parameters_name, &lbm, &net, &steer);
#else
  lbmReadParameters (input_parameters_name, &lbm, &net);
#endif


#ifdef STEER
  // create the derived datatype for the MPI_Bcast
  int steer_count = 25;
  int steer_blocklengths[25] = {1, 1, REG_MAX_NUM_STR_CMDS, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype steer_types[25] = {MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_DOUBLE,
				  MPI_DOUBLE, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER,
				  MPI_INTEGER, MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL,
				  MPI_REAL, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_REAL,
				  MPI_REAL, MPI_REAL, MPI_REAL, MPI_REAL};

  // calculate displacements
  MPI_Aint steer_disps[25];
  MPI_Datatype MPI_steer_type;
  
  
  steer_disps[0] = 0;
  
  for (int i = 1; i < steer_count; i++)
    {
      if (steer_types[i - 1] == MPI_INTEGER)
	{
	  steer_disps[i] = steer_disps[i - 1] + (sizeof(int) * steer_blocklengths[i - 1]);
	}
      else if (steer_types[i - 1] == MPI_DOUBLE)
	{
	  steer_disps[i] = steer_disps[i - 1] + (sizeof(double) * steer_blocklengths[i - 1]);
	}
      else if (steer_types[i - 1] == MPI_REAL)
	{
	  steer_disps[i] = steer_disps[i - 1] + (sizeof(float) * steer_blocklengths[i - 1]);
	}
    }
  
  MPI_Type_struct (steer_count, steer_blocklengths, steer_disps, steer_types, &MPI_steer_type);
  MPI_Type_commit (&MPI_steer_type);
  
  // initialize the steering library
  if(net.id == 0)
    {
      Steering_enable (REG_TRUE);
      
      reg_num_cmds = 2;
      reg_cmds[0] = REG_STR_STOP;
      reg_cmds[1] = REG_STR_PAUSE_INTERNAL;
      steer.status = Steering_initialize ("HemeLB", reg_num_cmds, reg_cmds);
    }

  // broadcast/collect status
  net.err = MPI_Bcast (&steer.status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

  // if broken, quit
  if(steer.status == REG_FAILURE)
    {
      net.err = MPI_Finalize ();
      return REG_FAILURE;
    }
#endif // STEER


#ifdef RG
  if(net.id == 0)
    {
      xdrSendBuffer = (char *)malloc(sizeof(char) * sizeToSend);
      
      pthread_mutex_init (&network_buffer_copy_lock, NULL);
      pthread_cond_init (&network_send_frame, NULL);
      
      pthread_attr_init (&pthread_attrib);
      pthread_attr_setdetachstate (&pthread_attrib, PTHREAD_CREATE_JOINABLE);
      
      pthread_create (&network_thread, &pthread_attrib, hemeLB_network, NULL);
      
      free (xdrSendBuffer);
    }
#endif // RG

  lbmInit (input_config_name, checkpoint_config_name, &lbm, &net);
  
  if (netFindTopology (&net) == 0)
    {
      fprintf (timings_ptr, "MPI_Attr_get failed, aborting\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  netInit (&lbm, &net, &rt);
  
#ifdef STEER
  rtReadParameters (rt_parameters_name, &net, &rt, &steer);
#else
  rtReadParameters (rt_parameters_name, &net, &rt);
#endif

  
  if (!lbm.is_checkpoint)
    {
      lbmSetOptimizedInitialConditions (&lbm, &net);
    }
  else
    {
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "Opening checkpoint file to read: %s\n", lbm.checkpoint_file_name);
	  fflush (timings_ptr);
	}
      void lbmSetInitialConditionsWithCheckpoint (LBM *lbm, Net *net);
    }
  
  rtInit (output_image_name, &net, &rt);
  
  stability = STABLE;
  checkpoint_count = 0;
  convergence_count = 0;
  ray_tracing_count = 0;
  
  
#ifdef STEER

  // register params with RealityGrid here
  if(net.id == 0)
    {
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
  if(steer.status == REG_FAILURE)
    {
      net.err = MPI_Finalize ();
      return REG_FAILURE;
    }

  reg_finished = 0;

  if(net.id == 0)
    {
      fprintf (timings_ptr, "STEER: RealityGrid library initialized and parameters registered.\n");
      fflush (timings_ptr);
    }
#endif // STEER
  
  net.err = MPI_Barrier (MPI_COMM_WORLD);
  
#ifndef BENCH
  simulation_time = myClock ();
  
  for (time_step = 1; time_step <= lbm.time_steps_max; time_step++)
    {
      write_checkpoint = 0;
      check_convergence = 0;
      perform_rt = 0;
      
      if (net.id == 0)
      	{
      	  fprintf (timings_ptr, "time step: %i, convergence error: %le\n",
		   time_step, fmin(1., lbm.convergence_error));
      	  fflush (timings_ptr);
      	}
      
#ifdef STEER
      // call steering control
      if (net.id == 0)
	{
	  steer.status = Steering_control (time_step,
					   &steer.num_params_changed,
					   steer_changed_param_labels,
					   &steer.num_recvd_cmds,
					   steer.recvd_cmds,
					   steer_recvd_cmd_params);
	}
      
      // broadcast/collect everything
      net.err = MPI_Bcast (&steer, 1, MPI_steer_type, 0, MPI_COMM_WORLD);

      if (steer.status != REG_SUCCESS)
	{
	  fprintf (timings_ptr, "STEER: I am %d and I detected that Steering_control failed.\n", net.id);
	  fflush (timings_ptr);
	  continue;
	}
      
      // process commands received
      for (int i = 0; i < steer.num_recvd_cmds; i++)
	{
	  switch (steer.recvd_cmds[i])
	    {
	    case REG_STR_STOP:
	      fprintf (timings_ptr, "STEER: I am %d and I've been told to STOP.\n", net.id);
	      fflush (timings_ptr);
	      reg_finished = 1;
	      break;
	    }
	} // end of command processing
      
      // process changed params
      // not bothered what changed, just copy across...
      
      if (steer.num_params_changed > 0)
	{
	  fprintf (timings_ptr, "STEER: I am %d and I was told that %d params changed.\n", net.id, steer.num_params_changed);
	  fflush (timings_ptr);
	  
	  lbmUpdateParameters (&lbm, &steer);
	  
	  rtUpdateParameters (&rt, &steer);
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
      
      // Between the rtRayTracingA/B calls, do not change any ray tracing
      // parameters.
      
      is_thread_locked = 0;
      
#ifdef RG
      if (net.id == 0 && perform_rt == 1)
	{
	  ///pthread_mutex_lock( &network_buffer_copy_lock ); ///
	  is_thread_locked = pthread_mutex_trylock ( &network_buffer_copy_lock );
	}
      net.err = MPI_Bcast (&is_thread_locked, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      
      if (perform_rt == 1 && is_thread_locked == 0)
      	{
      	  rtRayTracingA (AbsorptionCoefficients, &net, &rt);
      	}
      
      stability = lbmCycle (write_checkpoint, check_convergence, perform_rt,
			    &is_converged, &lbm, &net, &rt);
      
      if (perform_rt == 1 && is_thread_locked == 0)
	{
	  rtRayTracingB (AbsorptionCoefficients, &net, &rt);
	  
#ifdef RG
	  if (net.id == 0)
	    {
	      pthread_mutex_unlock (&network_buffer_copy_lock);
	      pthread_cond_signal (&network_send_frame);
	    }
#endif
	}
      
      if (stability == UNSTABLE || is_converged) break;

#ifdef STEER
      if (reg_finished == 1) break;
#endif // STEER
    }
  simulation_time = myClock () - simulation_time;
  
  time_step = min(time_step, lbm.time_steps_max);
#else // BENCH
  
  fluid_solver_time = myClock ();
  
  for (time_step = 1; time_step < 1000000000; time_step++)
    {
      if (net.id == 0)
	{
	  fprintf (timings_ptr, " fluid solver, time step: %i\n", time_step);
	  fflush (timings_ptr);
	}
      stability = lbmCycle (0, 0, 0, &is_converged, &lbm, &net, &rt);
      
      is_bench_section_finished = 0;
      
      if ((myClock () - fluid_solver_time) > minutes * (60. / 3.))
	{
	  is_bench_section_finished = 1;
	}
      net.err = MPI_Bcast (&is_bench_section_finished, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      
      if (is_bench_section_finished)
	{
	  break;
	}
    }
  fluid_solver_time = myClock () - fluid_solver_time;
  fluid_solver_time_steps = time_step;
  
  rt.image_frequency = 1;
  rt.flow_field_type = VELOCITY;
  rt.is_isosurface = 0;
  rt.cutoff = -EPSILON;
  fluid_solver_and_vr_time = myClock ();

  for (time_step = 1; time_step < 1000000000; time_step++)
    {
      if (net.id == 0)
	{
	  fprintf (timings_ptr, " fluid solver and vr, time step: %i\n", time_step);
	  fflush (timings_ptr);
	}
      rtRayTracingA (AbsorptionCoefficients, &net, &rt);
      
      stability = lbmCycle (0, 0, 1, &is_converged, &lbm, &net, &rt);
      
      rtRayTracingB (AbsorptionCoefficients, &net, &rt);
      
      is_bench_section_finished = 0;
      
      if ((myClock () - fluid_solver_and_vr_time) > minutes * (60. / 3.))
	{
	  is_bench_section_finished = 1;
	}
      net.err = MPI_Bcast (&is_bench_section_finished, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      
      if (is_bench_section_finished)
	{
	  break;
	}
    }
  fluid_solver_and_vr_time = myClock () - fluid_solver_and_vr_time;
  fluid_solver_and_is_time_steps = time_step;
  
  rt.image_frequency = 1;
  rt.flow_field_type = VELOCITY;
  rt.is_isosurface = 1;
  rt.cutoff = -EPSILON;
  fluid_solver_and_is_time = myClock ();
  
  for (time_step = 1; time_step < 1000000000; time_step++)
    {
     if (net.id == 0)
	{
	  fprintf (timings_ptr, " fluid solver and iso-surface, time step: %i\n", time_step);
	  fflush (timings_ptr);
	}
      rtRayTracingA (AbsorptionCoefficients, &net, &rt);
      
      stability = lbmCycle (0, 0, 1, &is_converged, &lbm, &net, &rt);
      
      rtRayTracingB (AbsorptionCoefficients, &net, &rt);
      
      is_bench_section_finished = 0;
      
      if ((myClock () - fluid_solver_and_is_time) > minutes * (60. / 3.))
	{
	  is_bench_section_finished = 1;
	}
      net.err = MPI_Bcast (&is_bench_section_finished, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      
      if (is_bench_section_finished)
	{
	  break;
	}
    }
  fluid_solver_and_is_time = myClock () - fluid_solver_and_is_time;
  fluid_solver_and_is_time_steps = time_step;
#endif // BENCH
  
#ifndef BENCH
  
  if (net.id == 0)
    {
      if (stability == STABLE)
  	{
  	  if (!is_converged)
  	    {
  	      fprintf (timings_ptr, " ATTENTION: SIMULATION NOT CONVERGED\n");
  	    }
  	}
      else
  	{
  	  fprintf (timings_ptr, " ATTENTION: INSTABILITY CONDITION OCCURRED\n");
  	  fprintf (timings_ptr, " AFTER %i time steps\n", time_step);
  	}
      fprintf (timings_ptr, "\n");
      fprintf (timings_ptr, "processors: %i, machines checked: %i\n", net.procs, net.machines);
      fprintf (timings_ptr, "time steps: %i \n", time_step);
      fprintf (timings_ptr, "time steps per second: %.3f\n\n", time_step / simulation_time);
    }
#else // BENCH
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "\n---------- BENCHMARKS RESULTS ----------\n");
      
      fprintf (timings_ptr, " procs checked: %i, machines checked: %i\n\n", net.procs, net.machines);
      fprintf (timings_ptr, " fluid sites: %i\n\n", lbm.total_fluid_sites);
      fprintf (timings_ptr, " time: %.3f, time steps per second: %.3f, MSUPS: %.3f\n\n",
	       fluid_solver_time, fluid_solver_time_steps / fluid_solver_time,
	       1.e-6 * lbm.total_fluid_sites / (fluid_solver_time / fluid_solver_time_steps));
      
      fprintf (timings_ptr, " time: %.3f, time steps per second with volume rendering: %.3f\n\n",
	       fluid_solver_and_vr_time, fluid_solver_and_vr_time_steps / fluid_solver_and_vr_time);
      
      fprintf (timings_ptr, " time: %.3f, time steps per second with isosurface: %.3f\n\n",
	       fluid_solver_and_is_time, fluid_solver_and_is_time_steps / fluid_solver_and_is_time);
    }
#endif
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "Opening output config file:\n %s\n\n", output_config_name);
      fflush (timings_ptr);
    }
  net.fo_time = myClock ();
  
  lbmWriteConfig (stability, output_config_name, 0, &lbm, &net);
  
  net.fo_time = myClock () - net.fo_time;
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "density  min, max: %le, %le\n", lbm.density_min, lbm.density_max);
      fprintf (timings_ptr, "velocity min, max: %le, %le\n", lbm.velocity_min, lbm.velocity_max);
      fprintf (timings_ptr, "stress   min, max: %le, %le\n", lbm.stress_min, lbm.stress_max);
      fprintf (timings_ptr, "\n");
      fprintf (timings_ptr, "domain decomposition time (s):             %.3f\n", net.dd_time);
      fprintf (timings_ptr, "pre-processing buffer management time (s): %.3f\n", net.bm_time);
      fprintf (timings_ptr, "input configuration reading time (s):      %.3f\n", net.fr_time);
      fprintf (timings_ptr, "flow field outputting time (s):            %.3f\n", net.fo_time);
      fflush (timings_ptr);
    }
  
  rtEnd (&net, &rt);
  lbmEnd (&lbm);
  netEnd (&net, &rt);

#ifdef RG
  if (net.id == 0)
    {
      // there are some problems if the following is called
      
      //pthread_join (network_thread, NULL);
    }
#endif // RG
  
  
#ifdef STEER
  if (net.id == 0)
    {
      Steering_finalize ();
      fprintf (timings_ptr, "STEER: Steering_finalize () called.\n");
      fflush (timings_ptr);
    }
#endif // STEER
  
  if (net.id == 0)
    {
      total_time = myClock () - total_time;
      fprintf (timings_ptr, "total time (s):                            %.3f\n", total_time);
      
      fclose (timings_ptr);
    }
  net.err = MPI_Finalize ();
  
  return(0);
}

