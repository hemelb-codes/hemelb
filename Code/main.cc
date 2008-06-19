// In this file, the functions useful to initiate/end the LB simulation
// and perform the dynamics are reported

#include "config.h"


#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#include <sys/stat.h>

#define MYPORT 65250
#define CONNECTION_BACKLOG 10


float steer_par[STEERABLE_PARAMETERS + 1];


FILE *timings_ptr;


void ColourPalette (float value, float col[])
{
  col[0] = fminf(1.F, value);
  col[1] = 0.;
  col[2] = fmaxf(0.F, 1.F - value);
}



char host_name[255];

// data per pixel inlude the data for the pixel location and 4 colours
//in RGB format, thus #bytes per pixel are (sizeof(int)+4
//rgb)=sizeof(int)+4*3*sizeof(unsigned char))
int bytes_per_pixel_data = sizeof(int) + 4 * sizeof(unsigned char);

// one int for colour_id and one for pixel id
u_int pixel_data_bytes = IMAGE_SIZE * bytes_per_pixel_data;

// it is assumed that the frame size is the only detail
u_int frame_details_bytes = 1 * sizeof(int);

char *xdrSendBuffer_pixel_data;
char *xdrSendBuffer_frame_details;



int recv_all (int sockid, char *buf, int *length)
{
  int received_bytes = 0;
  int bytes_left_to_receive = *length;
  int n;

  while (received_bytes < *length)
    {
      n = recv(sockid, buf+received_bytes, bytes_left_to_receive, 0);
      
      if (n == -1) break;
      
      received_bytes += n;
      bytes_left_to_receive -= n;
    }
  *length = received_bytes;
  
  return n == -1 ? -1 : 0;
}


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


void *hemeLB_steer (void *ptr)
{
  while(1) {

  long int read_fd = (long int)ptr;
  //printf("Kicking off steering thread with FD %i\n", (int)read_fd);
  
  int num_chars = STEERABLE_PARAMETERS * sizeof(float) / sizeof(char);
  int bytes = sizeof(char) * num_chars;
 
  char *xdr_steering_data = (char *)malloc(bytes);
  
  XDR xdr_steering_stream;
  
  
  xdrmem_create(&xdr_steering_stream, xdr_steering_data, bytes, XDR_DECODE);
 
  recv_all (read_fd, xdr_steering_data, &num_chars);
  
  for (int i = 0; i < STEERABLE_PARAMETERS; i++)
    {
      xdr_float (&xdr_steering_stream, &steer_par[i]);
    }
  free(xdr_steering_data);
  }
}


void *hemeLB_network (void *ptr)
{
  gethostname (host_name, 255);
  
  FILE *f = fopen ("env_details.asc","w");
  
  fprintf (f, "%s\n", host_name);
  fclose (f);
  
  fprintf (timings_ptr, "MPI 0 Hostname -> %s\n\n", host_name);
  
  
  int sock_fd;
  int new_fd;
  int yes = 1;
  
  int is_broken_pipe = 0;
  int frame_number = 0;
  
  pthread_t steering_thread;
  pthread_attr_t steering_thread_attrib; 
  
  
  signal(SIGPIPE, SIG_IGN); // Ignore a broken pipe 
  
  pthread_attr_init (&steering_thread_attrib);
  pthread_attr_setdetachstate (&steering_thread_attrib, PTHREAD_CREATE_JOINABLE);
  
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
      
      sin_size = sizeof (their_addr);
      
      if ((new_fd = accept (sock_fd, (struct sockaddr *)&their_addr, &sin_size)) == -1)
	{
	  perror("accept");
	  continue;
	}
      
      //fprintf (timings_ptr, "server: got connection from %s (FD %i)\n", inet_ntoa (their_addr.sin_addr), new_fd);
      printf ("RG thread: server: got connection from %s (FD %i)\n", inet_ntoa (their_addr.sin_addr), new_fd);
      
      pthread_create (&steering_thread, &steering_thread_attrib, hemeLB_steer, (void*)new_fd);	  
	  
      close(sock_fd);
      
      is_broken_pipe = 0;
      
      pthread_mutex_unlock ( &network_buffer_copy_lock );
      
      while (!is_broken_pipe)
	{
	  pthread_mutex_lock ( &network_buffer_copy_lock );
	  pthread_cond_wait (&network_send_frame, &network_buffer_copy_lock);
	  
	  int bytesSent = 0;
	  
	  XDR xdr_network_stream_frame_details;
	  XDR xdr_network_stream_pixel_data;
	  
	  
	  xdrmem_create (&xdr_network_stream_pixel_data, xdrSendBuffer_pixel_data,
			 pixel_data_bytes, XDR_ENCODE);
	  
	  xdrmem_create (&xdr_network_stream_frame_details, xdrSendBuffer_frame_details,
			 frame_details_bytes, XDR_ENCODE);
	  
	  for (int i = 0; i < col_pixels; i++)
	    {
	      xdrWritePixel (&col_pixel_recv[ i ], &xdr_network_stream_pixel_data, ColourPalette);
	    }
	  
	  int frameBytes = xdr_getpos(&xdr_network_stream_pixel_data);
	  
	  xdr_int (&xdr_network_stream_frame_details, &frameBytes);
	  
	  int detailsBytes = xdr_getpos(&xdr_network_stream_frame_details);
	  
	  int ret = send_all(new_fd, xdrSendBuffer_frame_details, &detailsBytes);
	  
          if (ret < 0) {
            is_broken_pipe = 1;
            break;
          } else {
            bytesSent += detailsBytes;
          }
	  
	  ret = send_all(new_fd, xdrSendBuffer_pixel_data, &frameBytes);
	  
          if (ret < 0) {
		    printf("RG thread: broken network pipe...\n");
            is_broken_pipe = 1;
			pthread_mutex_unlock ( &network_buffer_copy_lock );
            break;
          } else {
            bytesSent += frameBytes;
          }
	  
	  //fprintf (timings_ptr, "bytes sent %i\n", bytesSent);
	  printf ("RG thread: bytes sent %i\n", bytesSent);
	  
	  xdr_destroy (&xdr_network_stream_frame_details);
	  xdr_destroy (&xdr_network_stream_pixel_data);
	  
	  pthread_mutex_unlock ( &network_buffer_copy_lock );
	  
	  frame_number++;
	  
	} // while (is_broken_pipe == 0)
      
      close(new_fd);
      
    } // while(1)
}


void UpdateSteerableParameters (int *is_thread_locked, Vis *vis)
{
  steer_par[STEERABLE_PARAMETERS] = (float)(*is_thread_locked) + 0.1;
  
#ifndef NOMPI
  MPI_Bcast (steer_par, STEERABLE_PARAMETERS+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  *is_thread_locked = (int)steer_par[STEERABLE_PARAMETERS];
  
  if (*is_thread_locked) return;
  /*
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
  
  vis_velocity_threshold_max_inv = 1.F / velocity_max;
  vis_stress_threshold_max_inv   = 1.F / stress_max;
  */
}


int IsBenckSectionFinished (double minutes, double elapsed_time)
{
  int is_bench_section_finished = 0;
  
  
  if (elapsed_time > minutes * 60.)
    {
      is_bench_section_finished = 1;
    }
#ifndef NOMPI
  MPI_Bcast (&is_bench_section_finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  if (is_bench_section_finished)
    {
      return 1;
    }
  return 0;
}


void usage (char *progname)
{
  fprintf (timings_ptr, "Usage: %s path of the input files and minutes for benchmarking\n", progname);
  fprintf (timings_ptr, "if one wants to do a benchmark or\n");
  fprintf (timings_ptr, "number of pulsaticle cycles, time steps per cycle and\n");
  fprintf (timings_ptr, "voxel size in metres otherwise.\n");
  fprintf (timings_ptr, "The following files must be present in the path specified:\n");
  fprintf (timings_ptr, "config.dat, pars.asc rt_pars.asc\n");
}

int file_exists(const char * filename) {
        if (FILE * file = fopen(filename, "r")) {
                fclose(file);
                return 0;
    }
    return -1;
}

void check_file(const char * filename) {
        if(file_exists(filename) < 0 ) {
                fprintf(stderr,"Cannot open file %s\nExiting.\n", filename);
                exit(0);
        } else {
                fprintf(stderr,"Located file %s\n", filename);
	}
}

int main (int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output
  
  double simulation_time;
  double minutes;
  double fluid_solver_time;
  double fluid_solver_and_vis_time;
  double vis_without_compositing_time;
  
  int cycle_id;
  int total_time_steps, time_step, stability = STABLE;
  int depths;
  
  int fluid_solver_time_steps;
  int fluid_solver_and_vis_time_steps;
  int vis_without_compositing_time_steps;
  
  pthread_t network_thread;
  pthread_attr_t pthread_attrib;
  
  LBM lbm;
  
  Net net;
  
#ifndef NOMPI
  net.err = MPI_Init (&argc, &argv);
  net.err = MPI_Comm_size (MPI_COMM_WORLD, &net.procs);
  net.err = MPI_Comm_rank (MPI_COMM_WORLD, &net.id);
#else
  net.procs = 1;
  net.id = 0;
#endif
  
  if (argc == 3) // Check command line arguments
    {
      is_bench = 1;
      minutes = atof( argv[2] );
    }
  else if (argc == 5)
    {
      is_bench = 0;
      lbm.cycles_max = atoi( argv[2] );
      lbm.period     = atoi( argv[3] );
      lbm.voxel_size = atof( argv[4] );
    }
  else
    {
      if (net.id == 0) usage(argv[0]);

#ifndef NOMPI
      net.err = MPI_Abort (MPI_COMM_WORLD, 1);
      net.err = MPI_Finalize ();
#else
      exit(1);
#endif
    }

  double total_time = myClock();
  
  char* input_file_path( argv[1] );
  
  char input_config_name[256];
  char input_parameters_name[256];
  char output_config_name[256];
  char vis_parameters_name[256];
  char output_image_name[256];
  char timings_name[256];
  char procs_string[256];
  char image_name[256];
  char output_directory[256];
  
  strcpy ( input_config_name , input_file_path );
  strcat ( input_config_name , "/config.dat" );
  check_file(input_config_name);

  strcpy ( input_parameters_name , input_file_path );
  strcat ( input_parameters_name , "/pars.asc" );
  check_file(input_parameters_name);
  
  strcpy ( vis_parameters_name , input_file_path );
  strcat ( vis_parameters_name , "/rt_pars.asc" );
  check_file(vis_parameters_name);

  strcpy ( output_config_name , input_file_path );
  strcat ( output_config_name , "/out.dat" );
  
  /* Create directory for Images */
  strcpy(output_directory, input_file_path);
  strcat(output_directory, "/Images/");
  mkdir(output_directory, 0777);
  strcpy(output_directory, output_image_name);

  sprintf ( procs_string, "%i", net.procs);
  strcpy ( timings_name , input_file_path );
  strcat ( timings_name , "/timings" );
  strcat ( timings_name , procs_string );
  strcat ( timings_name , ".asc" );

  if (net.id == 0)
    {
      timings_ptr = fopen (timings_name, "w");
    }
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "***********************************************************\n");
      fprintf (timings_ptr, "Opening parameters file:\n %s\n", input_parameters_name);
      fprintf (timings_ptr, "Opening config file:\n %s\n", input_config_name);
      fprintf (timings_ptr, "Opening vis parameters file:\n %s\n\n", vis_parameters_name);
    }
  
  lbmReadParameters (input_parameters_name, &lbm, &net);
  
  if(net.id == 0)
    {
      xdrSendBuffer_pixel_data = (char *)malloc(pixel_data_bytes);
      xdrSendBuffer_frame_details = (char *)malloc(frame_details_bytes);
      
      pthread_mutex_init (&network_buffer_copy_lock, NULL);
      pthread_cond_init (&network_send_frame, NULL);
      
      pthread_attr_init (&pthread_attrib);
      pthread_attr_setdetachstate (&pthread_attrib, PTHREAD_CREATE_JOINABLE);
      
      pthread_create (&network_thread, &pthread_attrib, hemeLB_network, NULL);
    }
  
  lbmInit (input_config_name, &lbm, &net);
  
  if (netFindTopology (&net, &depths) == 0)
    {
      fprintf (timings_ptr, "MPI_Attr_get failed, aborting\n");
#ifndef NOMPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
  
  netInit (&lbm, &net);
  
  lbmSetInitialConditions (&net);
  
  visInit (&net, &vis);
  
  visReadParameters (vis_parameters_name, &lbm, &net, &vis);
  
  
  if (!is_bench)
    {
      int is_finished = 0;
      
      total_time_steps = 0;
      
      simulation_time = myClock ();
      
      for (cycle_id = 0; cycle_id < lbm.cycles_max && !is_finished; cycle_id++)
	{
	  for (time_step = 0; time_step < lbm.period; time_step++)
	    {
	      int perform_rt = 0;
	      int write_image = 0;
	      int stream_image = 0;
	      int is_thread_locked = 0;
	      
	      
	      ++total_time_steps;
	      
	      if ((time_step + 1)%vis_image_freq == 0)
		{
		  write_image = 1;
		}
	      if (total_time_steps%1 == 0)
		{
		  if (net.id == 0)
		    {
		      //pthread_mutex_lock( &network_buffer_copy_lock );
		      is_thread_locked = pthread_mutex_trylock ( &network_buffer_copy_lock );
		    }
		  UpdateSteerableParameters (&is_thread_locked, &vis);
		  
		  if (!is_thread_locked)
		    {
		      stream_image = 1;
		    }
		}
	      if (stream_image || write_image)
		{
		  perform_rt = 1;
		}
	      
	      // Between the visRenderA/B calls, do not change any vis
	      // parameters.
	      
	      if (perform_rt)
		{
		  visRenderA (ColourPalette, &net);
		}
	      lbmVaryBoundaryDensities (cycle_id, time_step, &lbm);
	      
	      stability = lbmCycle (cycle_id, time_step, perform_rt, &lbm, &net);
	      
	      if (write_image)
		{
		  char time_step_string[256];
		  
		  strcpy ( image_name , output_image_name ); /* At this point output_image_name is appended with Images */

		  int time_steps = time_step + 1;
		  
		  while (time_steps < 100000000)
		    {
		      strcat ( image_name , "0" );
		      time_steps *= 10;
		    } /* WTF? */
		  sprintf ( time_step_string, "%i", time_step + 1);
		  strcat ( image_name , time_step_string );
		  strcat ( image_name , ".dat" );
		}

	      if (perform_rt)
		{
		  visRenderB (write_image, image_name, ColourPalette, &net);
		}
	      if (net.id == 0)
		{
		  pthread_mutex_unlock (&network_buffer_copy_lock);
		  pthread_cond_signal (&network_send_frame);
		}
	      
	      if (stability == UNSTABLE)
		{
		  printf (" ATTENTION: INSTABILITY CONDITION OCCURRED\n");
		  printf (" AFTER %i total time steps\n", total_time_steps);
		  printf ("EXECUTION IS ABORTED\n");
#ifndef NOMPI
		  MPI_Abort (MPI_COMM_WORLD, 1);
#else
		  exit(1);
#endif
		  is_finished = 1;
		  break;
		}
	      if (net.id == 0)
		{
		  //printf ("time step: %i\n", time_step+1);
		}
	    }
	  if (net.id == 0)
	    {
	      //fprintf (timings_ptr, "cycle id: %i\n", cycle_id+1);
	      printf ("cycle id: %i\n", cycle_id+1);
	    }
	}
      simulation_time = myClock () - simulation_time;
      time_step = (1+min(time_step, lbm.period-1)) * min(cycle_id, lbm.cycles_max-1);
    }
  else // is_bench
    {
      double elapsed_time;
  
      int bench_period = (int)fmax(1., (1e+6 * net.procs) / lbm.total_fluid_sites);
      
      
      // benchmarking HemeLB's fluid solver only
      
      fluid_solver_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  stability = lbmCycle (0, 0, 0, &lbm, &net);
	  
	  // partial timings
	  elapsed_time = myClock () - fluid_solver_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " FS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      IsBenckSectionFinished (0.5, elapsed_time))
	    {
	      break;
	    }
	}
      fluid_solver_time_steps = (int)(time_step * minutes / (3. * 0.5) - time_step);
      fluid_solver_time = myClock ();
      
      for (time_step = 1; time_step <= fluid_solver_time_steps; time_step++)
	{
	  stability = lbmCycle (0, 0, 1, &lbm, &net);
	}
      fluid_solver_time = myClock () - fluid_solver_time;
      
      
      // benchmarking HemeLB's fluid solver and ray tracer
      
      vis_image_freq = 1;
      vis_compositing = 1;
      fluid_solver_and_vis_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  visRenderA (ColourPalette, &net);
	  
	  stability = lbmCycle (0, 0, 1, &lbm, &net);
	  
	  visRenderB (0, image_name, ColourPalette, &net);
	  
	  // partial timings
	  elapsed_time = myClock () - fluid_solver_and_vis_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " FS + VIS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      IsBenckSectionFinished (0.5, elapsed_time))
	    {
	      break;
	    }
	}
      fluid_solver_and_vis_time_steps = (int)(time_step * minutes / (3. * 0.5) - time_step);
      fluid_solver_and_vis_time = myClock ();
      
      for (time_step = 1; time_step <= fluid_solver_and_vis_time_steps; time_step++)
	{
	  visRenderA (ColourPalette, &net);
	  
	  stability = lbmCycle (0, 0, 1, &lbm, &net);
	  
	  visRenderB (0, image_name, ColourPalette, &net);
	}
      fluid_solver_and_vis_time = myClock () - fluid_solver_and_vis_time;
      
      
      
      // benchmarking HemeLB's ray tracer without compositing
      
      vis_compositing = 0;
      vis_without_compositing_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  visRenderA (ColourPalette, &net);
	  
	  // partial timings
	  elapsed_time = myClock () - vis_without_compositing_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " VIS - COMP, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      IsBenckSectionFinished (0.5, elapsed_time))
	    {
	      break;
	    }
	}
      vis_without_compositing_time_steps = (int)(time_step * minutes / (3. * 0.5) - time_step);
      vis_without_compositing_time = myClock ();
      
      for (time_step = 1; time_step <= vis_without_compositing_time_steps; time_step++)
	{
	  visRenderA (ColourPalette, &net);
	}
      vis_without_compositing_time = myClock () - vis_without_compositing_time;
    } // is_bench
  
  if (!is_bench)
    {  
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "\n");
	  fprintf (timings_ptr, "threads: %i, machines checked: %i\n\n", net.procs, net_machines);
	  fprintf (timings_ptr, "topology depths checked: %i\n\n", depths);
	  fprintf (timings_ptr, "fluid sites: %i\n\n", lbm.total_fluid_sites);
	  fprintf (timings_ptr, "cycles and total time steps: %i, %i \n\n", cycle_id, total_time_steps);
	  fprintf (timings_ptr, "time steps per second: %.3f\n\n", total_time_steps / simulation_time);
	}
    }
  else  // is_bench
    {
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "\n---------- BENCHMARK RESULTS ----------\n");
	  
	  fprintf (timings_ptr, "threads: %i, machines checked: %i\n\n", net.procs, net_machines);
	  fprintf (timings_ptr, "topology depths checked: %i\n\n", depths);
	  fprintf (timings_ptr, "fluid sites: %i\n\n", lbm.total_fluid_sites);
	  fprintf (timings_ptr, " FS, time steps per second: %.3f, MSUPS: %.3f, time: %.3f\n\n",
		   fluid_solver_time_steps / fluid_solver_time,
		   1.e-6 * lbm.total_fluid_sites / (fluid_solver_time / fluid_solver_time_steps),
		   fluid_solver_time);
	  
	  fprintf (timings_ptr, " FS + VIS, time steps per second: %.3f, time: %.3f\n\n",
		   fluid_solver_and_vis_time_steps / fluid_solver_and_vis_time, fluid_solver_and_vis_time);
	  
	  fprintf (timings_ptr, " VR - COMP, time steps per second: %.3f, time: %.3f\n\n",
		   vis_without_compositing_time_steps / vis_without_compositing_time, vis_without_compositing_time);
	}
    }
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "Opening output config file:\n %s\n\n", output_config_name);
      fflush (timings_ptr);
    }
  net.fo_time = myClock ();
  
  lbmWriteConfig (stability, output_config_name, &lbm, &net);
  
  net.fo_time = myClock () - net.fo_time;
  
  if (net.id == 0)
    {
      if (!is_bench)
	{
	  double pressure_min = lbmConvertPressureToPhysicalUnits (lbm_density_min * Cs2, &lbm);
	  double pressure_max = lbmConvertPressureToPhysicalUnits (lbm_density_max * Cs2, &lbm);
	  
	  double velocity_min = lbmConvertVelocityToPhysicalUnits (lbm_velocity_min, &lbm);
	  double velocity_max = lbmConvertVelocityToPhysicalUnits (lbm_velocity_max, &lbm);
	  
	  double stress_min = lbmConvertStressToPhysicalUnits (lbm_stress_min, &lbm);
	  double stress_max = lbmConvertStressToPhysicalUnits (lbm_stress_max, &lbm);
	  
	  fprintf (timings_ptr, "pressure min, max (mmHg): %le, %le\n", pressure_min, pressure_max);
	  fprintf (timings_ptr, "velocity min, max (m/s) : %le, %le\n", velocity_min, velocity_max);
	  fprintf (timings_ptr, "stress   min, max (Pa)  : %le, %le\n", stress_min, stress_max);
	}
      fprintf (timings_ptr, "\n");
      fprintf (timings_ptr, "domain decomposition time (s):             %.3f\n", net.dd_time);
      fprintf (timings_ptr, "pre-processing buffer management time (s): %.3f\n", net.bm_time);
      fprintf (timings_ptr, "input configuration reading time (s):      %.3f\n", net.fr_time);
      fprintf (timings_ptr, "flow field outputting time (s):            %.3f\n", net.fo_time);
      
      total_time = myClock () - total_time;
      fprintf (timings_ptr, "total time (s):                            %.3f\n\n", total_time);
      
      fprintf (timings_ptr, "Sub-domains info:\n\n");
      
      for (int n = 0; n < net.procs; n++)
	{
	  fprintf (timings_ptr, "rank: %i, fluid sites: %i\n", n, net.fluid_sites[ n ]);
	}
      
      fclose (timings_ptr);
    }
  
  visEnd ();
  netEnd (&net);
  lbmEnd (&lbm);
  
  if (net.id == 0)
    {
      // there are some problems if the following function is called
      
      //pthread_join (network_thread, NULL);
      free(xdrSendBuffer_frame_details);
      free(xdrSendBuffer_pixel_data);
    }
  
  
#ifndef NOMPI
  net.err = MPI_Finalize ();
#endif
  
  return(0);
}
