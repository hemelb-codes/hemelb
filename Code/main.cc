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

float globalLongitude = 0.;


void visUpdateLongitude (char *parameters_file_name, LBM *lbm, Net *net, Vis *vis)
{
  FILE *parameters_file;

  float par_to_send[14];
  float ctr_x, ctr_y, ctr_z;
  float longitude, latitude;
  float zoom;
  float density_max, velocity_max, stress_max;
  float dummy;
  
  
  if (net->id == 0)
    {
      parameters_file = fopen (parameters_file_name, "r");
      
      fscanf (parameters_file, "%e \n", &dummy);
      fscanf (parameters_file, "%e \n", &dummy);
      fscanf (parameters_file, "%e \n", &ctr_x);
      fscanf (parameters_file, "%e \n", &ctr_y);
      fscanf (parameters_file, "%e \n", &ctr_z);
      fscanf (parameters_file, "%e \n", &longitude);
      fscanf (parameters_file, "%e \n", &latitude);
      fscanf (parameters_file, "%e \n", &zoom);
      
      fscanf (parameters_file, "%i \n", &vis->image_freq);
      fscanf (parameters_file, "%i \n", &lbm->flow_field_type);
      fscanf (parameters_file, "%i \n", &vis->mode);
      fscanf (parameters_file, "%e \n", &vis->absorption_factor);
      fscanf (parameters_file, "%e \n", &vis->cutoff);
      fscanf (parameters_file, "%e \n", &density_max);
      fscanf (parameters_file, "%e \n", &velocity_max);
      fscanf (parameters_file, "%e \n", &stress_max);

      fclose (parameters_file);
      
      par_to_send[  0 ] = ctr_x;
      par_to_send[  1 ] = ctr_y;
      par_to_send[  2 ] = ctr_z;
      par_to_send[  3 ] = longitude;
      par_to_send[  4 ] = latitude;
      par_to_send[  5 ] = zoom;
      par_to_send[  6 ] = 0.1 + (float)vis->image_freq;
      par_to_send[  7 ] = 0.1 + (float)lbm->flow_field_type;
      par_to_send[  8 ] = 0.1 + (float)vis->mode;
      par_to_send[  9 ] = vis->absorption_factor;
      par_to_send[ 10 ] = vis->cutoff;
      par_to_send[ 11 ] = density_max;
      par_to_send[ 12 ] = velocity_max;
      par_to_send[ 13 ] = stress_max;
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 14, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  ctr_x                  =      par_to_send[  0 ];
  ctr_y                  =      par_to_send[  1 ];
  ctr_z                  =      par_to_send[  2 ];
  longitude              =      par_to_send[  3 ];
  latitude               =      par_to_send[  4 ];
  zoom                   =      par_to_send[  5 ];
  vis->image_freq        = (int)par_to_send[  6 ];
  lbm->flow_field_type   = (int)par_to_send[  7 ];
  vis->mode              = (int)par_to_send[  8 ];
  vis->absorption_factor =      par_to_send[  9 ];
  vis->cutoff            =      par_to_send[ 10 ];
  density_max            =      par_to_send[ 11 ];
  velocity_max           =      par_to_send[ 12 ];
  stress_max             =      par_to_send[ 13 ];
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
		 PIXELS_X, PIXELS_Y,
		 ctr_x, ctr_y, ctr_z,
		 5.F * vis->system_size,
		 globalLongitude, latitude,
		 0.5F * (5.F * vis->system_size),
		 zoom);
  
  if (lbm->flow_field_type == DENSITY)
    {
      vis->flow_field_value_max_inv = 1.F / density_max;
    }
  else if (lbm->flow_field_type == VELOCITY)
    {
      vis->flow_field_value_max_inv = 1.F / velocity_max;
    }
  else
    {
      vis->flow_field_value_max_inv = 1.F / stress_max;
    }
  vis->velocity_max_inv = 1.F / velocity_max;
}



#ifdef RG

char host_name[255];

// data per pixel are colour id and pixel id (2 * sizeof(int) bytes)
int data_per_pixel = 2;
int bytes_per_pixel_data = data_per_pixel * sizeof(int);

// one int for colour_id and one for pixel id
u_int pixel_data_bytes = IMAGE_SIZE * bytes_per_pixel_data;

// it is assumed that the frame size is the only detail
u_int frame_details_bytes = 1 * sizeof(int);

char *xdrSendBuffer_pixel_data;
char *xdrSendBuffer_frame_details;

int bits_per_char = sizeof(char) * 8;
int bits_per_two_chars = 2 * bits_per_char;


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

// #ifndef STEER
  FILE *f = fopen ("env_details.asc","w");
  
  fprintf (f, "%s\n", host_name);
  fclose (f);
  
  fprintf (timings_ptr, "MPI 0 Hostname -> %s\n\n", host_name);
// #endif
  
  signal(SIGPIPE, SIG_IGN); // Ignore a broken pipe
  
  int sock_fd;
  int new_fd;
  int yes = 1;
  
  int is_broken_pipe = 0;
  int frame_number = 0;
  
  int pixel_r, pixel_g, pixel_b;
  int pixel_i, pixel_j;
  int colour_id, pixel_id;
  
  ColPixel *col_pixel_p;
  
  
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
      printf ("server: got connection from %s\n", inet_ntoa (their_addr.sin_addr));
      
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
	  
	  for (int i = 0; i < vis.col_pixels_locked; i++)
	    {
	      col_pixel_p = &vis.col_pixel_locked[ i ];
	      
	      pixel_r = max(0, min(255, (int)(255.F - col_pixel_p->r)));
	      pixel_g = max(0, min(255, (int)(255.F - col_pixel_p->g)));
	      pixel_b = max(0, min(255, (int)(255.F - col_pixel_p->b)));
	      
	      pixel_i = col_pixel_p->i;
	      pixel_j = col_pixel_p->j;
	      
	      colour_id = (pixel_r << bits_per_two_chars) + (pixel_g << bits_per_char) + pixel_b;
	      pixel_id = (pixel_i << bits_per_two_chars) + pixel_j;
	      
	      xdr_int (&xdr_network_stream_pixel_data, &colour_id);
	      xdr_int (&xdr_network_stream_pixel_data, &pixel_id);
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
            is_broken_pipe = 1;
            break;
          } else {
            bytesSent += frameBytes;
          }
	  
	  fprintf (timings_ptr, "bytes sent %i\n", bytesSent);
	  printf ("bytes sent %i\n", bytesSent);
	  
	  xdr_destroy (&xdr_network_stream_frame_details);
	  xdr_destroy (&xdr_network_stream_pixel_data);
	  
	  pthread_mutex_unlock ( &network_buffer_copy_lock );
	  
	  frame_number++;
	  
	} // while (is_broken_pipe == 0)
      
      close(new_fd);
      
    } // while(1)
}

#endif // RG


inline void ColourPalette (float value, float col[])
{
  value = fminf(1.F, value);
  
  col[0] = (1.F - value) * (1.F - value);
  col[1] = value * (1.F - value);
  col[2] = value * value;
}


void TimeVaryingDensities (int period, int time_step, int inlets, int outlets,
			   double inlet_density[], double outlet_density[])
{
  double K = 0.00025;
  double w = 2. * PI / period;
  
  inlet_density[0]  = 1. + (((32.*0.5/Cs2) * K) * cos(w * (double)time_step + 0.5 * PI));
  outlet_density[0] = 1. - (((32.*0.5/Cs2) * K) * cos(w * (double)time_step + 0.5 * PI));
}


int IsBenckSectionFinished (double minutes, double elapsed_time)
{
  int is_bench_section_finished = 0;
#ifndef NOMPI
  int err;
#endif
  
  if (elapsed_time > minutes * 60.)
    {
      is_bench_section_finished = 1;
    }
#ifndef NOMPI
  err = MPI_Bcast (&is_bench_section_finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  if (is_bench_section_finished)
    {
      return 1;
    }
  return 0;
}


void usage (char *progname)
{
#ifndef BENCH
  fprintf (timings_ptr, "Usage: %s path of the input files\n", progname);
  fprintf (timings_ptr, "the following files must be present: config.dat, pars.asc\n");
  fprintf (timings_ptr, "check.dat (if it is a checkpoint file and must be used), vis_pars.asc\n");
#else
  fprintf (timings_ptr, "Usage: %s path of the input files and minutes for benchmarking\n", progname);
  fprintf (timings_ptr, "the following files must be present: config.dat, pars.asc\n");
  fprintf (timings_ptr, "check.dat (if it is a checkpoint file and must be used), vis_pars.asc\n");
#endif
}


int main (int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output
  
#ifndef BENCH
  double simulation_time;
#else
  double minutes;
  double fluid_solver_time, fluid_solver_and_vr_time, fluid_solver_and_is_time;
#endif // BENCH
  
#ifndef BENCH
  int cycle_id;
#endif
  int time_step, stability, is_converged;
  int *proc_fluid_sites;
  int depths;
  
#ifdef BENCH
  int fluid_solver_time_steps;
  int fluid_solver_and_vr_time_steps;
  int fluid_solver_and_is_time_steps;
#else
#ifndef TD
  int checkpoint_count = 0;
#endif
  int conv_count = 0;
  int ray_tracing_count = 0;
#ifndef TD
  int write_checkpoint;
#endif
  int check_conv, perform_vis;
  int is_thread_locked;
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
  
  
#ifndef NOMPI
  net.err = MPI_Init (&argc, &argv);
  net.err = MPI_Comm_size (MPI_COMM_WORLD, &net.procs);
  net.err = MPI_Comm_rank (MPI_COMM_WORLD, &net.id);
#else
  net.procs = 1;
  net.id = 0;
#endif
  
  double total_time = myClock ();
  
  char *input_file_path( argv[1] );
  
  char input_config_name[256];
  char input_parameters_name[256];
  char output_config_name[256];
  char checkpoint_config_name[256];
  char vis_parameters_name[256];
  char output_image_name[256];
  char timings_name[256];
  
  
  strcpy ( input_config_name , input_file_path );
  strcat ( input_config_name , "/config.dat" );

  strcpy ( input_parameters_name , input_file_path );
  strcat ( input_parameters_name , "/pars.asc" );
  
  strcpy ( output_config_name , input_file_path );
  strcat ( output_config_name , "/out.dat" );
  
  strcpy ( checkpoint_config_name , input_file_path );
  strcat ( checkpoint_config_name , "/check.dat" );
  
  strcpy ( vis_parameters_name , input_file_path );
  strcat ( vis_parameters_name , "/rt_pars.asc" );
  
  strcpy ( output_image_name , input_file_path );
  strcat ( output_image_name , "/image.dat" );
  
  strcpy ( timings_name , input_file_path );
  strcat ( timings_name , "/timings.asc" );
  
  if (net.id == 0)
    {
      timings_ptr = fopen (timings_name, "w");
    }
  
#ifdef BENCH
  if (argc != 3)
#else
  if (argc != 2 && argc != 3)
#endif
    {
      if (net.id == 0) usage(argv[0]);
      
#ifndef NOMPI
      net.err = MPI_Abort (MPI_COMM_WORLD, 1);
      net.err = MPI_Finalize ();
#else
      exit(1);
#endif
    }
  
#ifdef BENCH
  minutes = atof( argv[2] );
#endif
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "***********************************************************\n");
      fprintf (timings_ptr, "Opening parameters file:\n %s\n", input_parameters_name);
      fprintf (timings_ptr, "Opening config file:\n %s\n", input_config_name);
      fprintf (timings_ptr, "Opening vis parameters file:\n %s\n\n", vis_parameters_name);
    }
  
  lbm.inlet_density = NULL;
  lbm.outlet_density = NULL;
  
#ifdef STEER
  lbmReadParameters (input_parameters_name, &lbm, &net, &steer);
#else
  lbmReadParameters (input_parameters_name, &lbm, &net);
#endif
  
  
#ifdef RG
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
#endif // RG
  
  
#ifdef STEER
  // create the derived datatype for the MPI_Bcast
  int steer_count = 24;
  int steer_blocklengths[24] = {1, 1, REG_MAX_NUM_STR_CMDS, 1,
				1, 1, 1, 1, 1,
				1, 1, 1,
				1, 1, 1,
				1, 1, 1,
				1, 1, 1, 1, 1, 1, 1,
				1};
#ifndef NOMPI
  MPI_Datatype steer_types[24] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT,
				  MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT,
				  MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				  MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				  MPI_INT, MPI_INT, MPI_INT,
				  MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				  MPI_UB};
  
  MPI_Aint steer_disps[24];
  MPI_Datatype MPI_steer_type;
  
  // calculate displacements
  
  steer_disps[0] = 0;
  
  for(int i = 1; i < steer_count; i++) {
    switch(steer_types[i - 1]) {
    case MPI_INT:
      steer_disps[i] = steer_disps[i - 1] + (sizeof(int) * steer_blocklengths[i - 1]);
      break;
    case MPI_DOUBLE:
      steer_disps[i] = steer_disps[i - 1] + (sizeof(double) * steer_blocklengths[i - 1]);
      break;
    case MPI_FLOAT:
      steer_disps[i] = steer_disps[i - 1] + (sizeof(float) * steer_blocklengths[i - 1]);
      break;
    }
  }
  
  MPI_Type_struct (steer_count, steer_blocklengths, steer_disps, steer_types, &MPI_steer_type);
  MPI_Type_commit (&MPI_steer_type);
#endif
  
  
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
#ifndef NOMPI
  net.err = MPI_Bcast (&steer.status, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  // if broken, quit
  if(steer.status == REG_FAILURE)
    {
#ifndef NOMPI
      net.err = MPI_Finalize ();
#endif
      return REG_FAILURE;
    }
#endif // STEER
  
  proc_fluid_sites = (int *)malloc(sizeof(int) * net.procs);
  
  lbmInit (input_config_name, checkpoint_config_name, &lbm, &net);
  
  if (netFindTopology (&net, &depths) == 0)
    {
      fprintf (timings_ptr, "MPI_Attr_get failed, aborting\n");
#ifndef NOMPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
  
  netInit (&lbm, &net, proc_fluid_sites);
  
  if (!lbm.is_checkpoint)
    {
      lbmSetInitialConditions (&lbm, &net);
    }
  else
    {
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "Opening checkpoint file to read: %s\n", lbm.checkpoint_file_name);
	  fflush (timings_ptr);
	}
      lbmSetInitialConditionsWithCheckpoint (&lbm, &net);
    }
  
  visInit (output_image_name, &net, &vis);
  
  stability = STABLE;
  
#ifdef STEER

  // register params with RealityGrid here
  if(net.id == 0)
    {
      // read only and only if displaying
//#ifdef RG
//      steer.status = Register_param("Display Host", REG_FALSE,
//				    (void*)(&host_name), REG_CHAR, "", "");
//#endif 
      
      // LBM params
      steer.status = Register_param("Tau", REG_TRUE,
				    (void*)(&steer.tau), REG_DBL, "0.5", "");
      steer.status = Register_param("Tolerance", REG_TRUE,
				    (void*)(&steer.tolerance), REG_DBL, "0.0", "0.1");
      steer.status = Register_param("Max time steps",REG_TRUE,
				    (void*)(&steer.max_cycles), REG_INT, "1", "");
      steer.status = Register_param("Conv frequency", REG_TRUE,
				    (void*)(&steer.conv_freq), REG_INT, "1", "");
#ifndef TD
      steer.status = Register_param("Checkpoint frequency", REG_TRUE,
				    (void*)(&steer.check_freq), REG_INT, "1", "");
#else
      steer.status = Register_param("Checkpoint frequency", REG_TRUE,
				    (void*)(&steer.period), REG_INT, "1", "");
#endif
      // Vis params
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
      steer.status = Register_param("Vis mode", REG_TRUE,
				    (void*)(&steer.mode), REG_INT, "0", "2");
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
  
#ifndef NOMPI
  // broadcast/collect status
  net.err = MPI_Bcast(&steer.status, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  // if broken, quit
  if(steer.status == REG_FAILURE)
    {
#ifndef NOMPI
      net.err = MPI_Finalize ();
#endif
      return REG_FAILURE;
    }

  reg_finished = 0;

  if(net.id == 0)
    {
      fprintf (timings_ptr, "STEER: RealityGrid library initialized and parameters registered.\n");
      fflush (timings_ptr);
    }
#endif // STEER
  
#ifdef STEER
  visReadParameters (vis_parameters_name, &lbm, &net, &vis, &steer);
#else
  visReadParameters (vis_parameters_name, &lbm, &net, &vis);
#endif
  
#ifdef STEER
  int steering_time_step = 0;
#endif
  
#ifndef BENCH
  int is_stop = 0;
  
  simulation_time = myClock ();
  
  for (cycle_id = 0; cycle_id < lbm.cycles_max && !is_stop; cycle_id++)
    {
      for (time_step = 0; time_step < lbm.period; time_step++)
	{
	  //globalLongitude += 1.F;
	  //visUpdateLongitude (vis_parameters_name, &lbm, &net, &vis);
	  
#ifndef TD
	  write_checkpoint = 0;
#endif
	  check_conv = 0;
	  perform_vis = 0;
	  
#ifdef STEER
	  // call steering control
	  if (net.id == 0)
	    {
	      steer.status = Steering_control (++steering_time_step,
					       &steer.num_params_changed,
					       steer_changed_param_labels,
					       &steer.num_recvd_cmds,
					       steer.recvd_cmds,
					       steer_recvd_cmd_params);
	    }
	  
#ifndef NOMPI
	  // broadcast/collect everything
	  net.err = MPI_Bcast (&steer, 1, MPI_steer_type, 0, MPI_COMM_WORLD);
#endif
	  if (steer.status != REG_SUCCESS)
	    {
	      fprintf (stderr, "STEER: I am %d and I detected that Steering_control failed.\n", net.id);
	      fflush(stderr);  
	      continue;
	    }
	  
	  // process commands received
	  for (int i = 0; i < steer.num_recvd_cmds; i++)
	    {
	      switch (steer.recvd_cmds[i])
		{
		case REG_STR_STOP:
		  fprintf (stderr, "STEER: I am %d and I've been told to STOP.\n", net.id);
		  fflush(stderr);
		  reg_finished = 1;
		  break;
		}
	    } // end of command processing
	  
	  // process changed params
	  // not bothered what changed, just copy across...
	  
	  if (steer.num_params_changed > 0)
	    {
	      fprintf (stderr, "STEER: I am %d and I was told that %d params changed.\n", net.id, steer.num_params_changed);
	      fflush(stderr);
	      
	      lbmUpdateParameters (&lbm, &steer);
	      
	      visUpdateParameters (&vis, &steer);
	    }
	  // end of param processing
	  
#endif // STEER
	  
#ifndef TD
	  if (++checkpoint_count >= lbm.checkpoint_freq)
	    {
	      write_checkpoint = 1;
	      checkpoint_count = 0;
	    }
#endif
	  if (++conv_count >= lbm.conv_freq)
	    {
	      check_conv = 1;
	      conv_count = 0;
	    }
	  if (++ray_tracing_count >= vis.image_freq)
	    {
	      perform_vis = 1;
	      ray_tracing_count = 0;
	    }
	  
	  // Between the visRenderA/B calls, do not change any vis
	  // parameters.
	  
	  is_thread_locked = 0;
	  
#ifdef RG
	  if (net.id == 0 && perform_vis == 1)
	    {
	      ///pthread_mutex_lock( &network_buffer_copy_lock ); ///
	      is_thread_locked = pthread_mutex_trylock ( &network_buffer_copy_lock );
	    }
#ifndef NOMPI
	  net.err = MPI_Bcast (&is_thread_locked, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#endif // RG
	  
	  if (perform_vis == 1 && is_thread_locked == 0)
	    {
	      visRenderA (ColourPalette, &net, &vis);
	    }
#ifndef TD
	  stability = lbmCycle (write_checkpoint, check_conv, perform_vis, &is_converged, &lbm, &net);
#else
	  TimeVaryingDensities (lbm.period, time_step, lbm.inlets, lbm.outlets,
	  			lbm.inlet_density, lbm.outlet_density);
	  
	  stability = lbmCycle (cycle_id, time_step, check_conv, perform_vis, &is_converged, &lbm, &net);
#endif // TD
	  if (perform_vis == 1 && is_thread_locked == 0)
	    {
	      visRenderB (&net, &vis);
#ifdef RG
	      if (net.id == 0)
		{
		  pthread_mutex_unlock (&network_buffer_copy_lock);
		  pthread_cond_signal (&network_send_frame);
		}
#endif
	    }
#ifdef TD
	  if (net.id == 0)
	    {
	      fprintf (timings_ptr, "time step: %i\n", time_step+1);
	      printf ("time step: %i\n", time_step+1);
	    }
#endif
	  if (stability == UNSTABLE || is_converged)
	    {
	      is_stop = 1;
	      break;
	    }
#ifdef STEER
	  if (reg_finished == 1)
	    {
	      is_stop = 1;
	      break;
	    }
#endif // STEER
	}
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "cycle id: %i, conv error: %le\n",
		   cycle_id+1, fmin(1., lbm.conv_error));
	  printf ("cycle id: %i, conv error: %le\n",
		  cycle_id+1, fmin(1., lbm.conv_error));
	}
    }
  simulation_time = myClock () - simulation_time;
  time_step = (1+min(time_step, lbm.period-1)) * min(cycle_id, lbm.cycles_max-1);
#else // BENCH
  
  double elapsed_time;
  
  // benchmarking HemeLB's fluid solver only
  
  fluid_solver_time = myClock ();
  
  for (time_step = 1; time_step <= 1000000000; time_step++)
    {
#ifndef TD
      stability = lbmCycle (0, 0, 0, &is_converged, &lbm, &net);
#else
      stability = lbmCycle (0, 0, 0, 0, &is_converged, &lbm, &net);
#endif // TD
      
      // partial timings
      elapsed_time = myClock () - fluid_solver_time;
      
      if (time_step%100 == 1 && net.id == 0)
	{
	  fprintf (timings_ptr, " FS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		   elapsed_time, time_step, time_step / elapsed_time);
	  fprintf (stderr, " FS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		   elapsed_time, time_step, time_step / elapsed_time);
	}
      if (time_step%100 == 1 &&
	  IsBenckSectionFinished (minutes / 3., elapsed_time))
	{
	  break;
	}
    }
  fluid_solver_time = myClock () - fluid_solver_time;
  fluid_solver_time_steps = time_step;
  
  
  // benchmarking HemeLB's fluid solver and volume rendering
  
  lbm.flow_field_type = VELOCITY;
  vis.image_freq = 1;
  vis.mode = 0;
  vis.cutoff = -EPSILON;
  fluid_solver_and_vr_time = myClock ();
  
  for (time_step = 1; time_step <= 1000000000; time_step++)
    {
      visRenderA (ColourPalette, &net, &vis);
#ifndef TD
      stability = lbmCycle (0, 0, 0, &is_converged, &lbm, &net);
#else
      stability = lbmCycle (0, 0, 0, 0, &is_converged, &lbm, &net);
#endif // TD
      visRenderB (&net, &vis);
      
      // partial timings
      elapsed_time = myClock () - fluid_solver_and_vr_time;
      
      if (time_step%100 == 1 && net.id == 0)
	{
	  fprintf (timings_ptr, " FS + VR, time: %.3f, time step: %i, time steps/s: %.3f\n",
		   elapsed_time, time_step, time_step / elapsed_time);
	  fprintf (stderr, " FS + VR, time: %.3f, time step: %i, time steps/s: %.3f\n",
		   elapsed_time, time_step, time_step / elapsed_time);
	}
      if (time_step%100 == 1 &&
	  IsBenckSectionFinished (minutes / 3., elapsed_time))
	{
	  break;
	}
    }
  fluid_solver_and_vr_time = myClock () - fluid_solver_and_vr_time;
  fluid_solver_and_vr_time_steps = time_step;
  
  
  // benchmarking HemeLB's fluid solver and iso-surface
  
  lbm.flow_field_type = VELOCITY;
  vis.image_freq = 1;
  vis.mode = 1;
  vis.cutoff = -EPSILON;
  fluid_solver_and_is_time = myClock ();
  
  for (time_step = 1; time_step <= 1000000000; time_step++)
    {
      visRenderA (ColourPalette, &net, &vis);
#ifndef TD
      stability = lbmCycle (0, 0, 0, &is_converged, &lbm, &net);
#else
      stability = lbmCycle (0, 0, 0, 0, &is_converged, &lbm, &net);
#endif // TD
      visRenderB (&net, &vis);
      
      // partial timings
      elapsed_time = myClock () - fluid_solver_and_is_time;
      
      if (time_step%100 == 1 && net.id == 0)
	{
	  fprintf (timings_ptr, " FS + IS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		   elapsed_time, time_step, time_step / elapsed_time);
	  fprintf (stderr, " FS + IS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		   elapsed_time, time_step, time_step / elapsed_time);
	}
      if (time_step%100 == 1 &&
	  IsBenckSectionFinished (minutes / 3., elapsed_time))
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
      fprintf (timings_ptr, "processors: %i, machines checked: %i\n\n", net.procs, net.machines);
      fprintf (timings_ptr, "topology depths checked: %i\n\n", depths);
      fprintf (timings_ptr, "fluid sites: %i\n\n", lbm.total_fluid_sites);
      fprintf (timings_ptr, "time steps: %i \n\n", time_step);
      fprintf (timings_ptr, "time steps per second: %.3f\n\n", time_step / simulation_time);
    }
#else // BENCH
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "\n---------- BENCHMARK RESULTS ----------\n");
      
      fprintf (timings_ptr, "procs checked: %i, machines checked: %i\n\n", net.procs, net.machines);
      fprintf (timings_ptr, "topology depths checked: %i\n\n", depths);
      fprintf (timings_ptr, "fluid sites: %i\n\n", lbm.total_fluid_sites);
      fprintf (timings_ptr, "time steps: %i \n\n", time_step);
      fprintf (timings_ptr, "time steps per second: %.3f, MSUPS: %.3f, time: %.3f\n\n",
	       fluid_solver_time_steps / fluid_solver_time,
	       1.e-6 * lbm.total_fluid_sites / (fluid_solver_time / fluid_solver_time_steps),
	       fluid_solver_time);
      
      fprintf (timings_ptr, "time steps per second with volume rendering: %.3f, time: %.3f\n\n",
	       fluid_solver_and_vr_time_steps / fluid_solver_and_vr_time, fluid_solver_and_vr_time);
      
      fprintf (timings_ptr, "time steps per second with isosurface: %.3f, time: %.3f\n\n",
	       fluid_solver_and_is_time_steps / fluid_solver_and_is_time, fluid_solver_and_is_time);
    }
#endif
  
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
      fprintf (timings_ptr, "density  min, max: %le, %le\n", lbm.density_min, lbm.density_max);
      fprintf (timings_ptr, "velocity min, max: %le, %le\n", lbm.velocity_min, lbm.velocity_max);
      fprintf (timings_ptr, "stress   min, max: %le, %le\n", lbm.stress_min, lbm.stress_max);
      fprintf (timings_ptr, "\n");
      fprintf (timings_ptr, "domain decomposition time (s):             %.3f\n", net.dd_time);
      fprintf (timings_ptr, "pre-processing buffer management time (s): %.3f\n", net.bm_time);
      fprintf (timings_ptr, "input configuration reading time (s):      %.3f\n", net.fr_time);
      fprintf (timings_ptr, "flow field outputting time (s):            %.3f\n", net.fo_time);
    }
  
  visEnd (&net, &vis);
  netEnd (&net);
  lbmEnd (&lbm);
  
#ifdef RG
  if (net.id == 0)
    {
      // there are some problems if the following function is called
      
      //pthread_join (network_thread, NULL);
      free(xdrSendBuffer_frame_details);
      free(xdrSendBuffer_pixel_data);
    }
#endif // RG
  
  
#ifdef STEER
  if (net.id == 0)
    {
      Steering_finalize ();
      
      fprintf (stderr, "STEER: Steering_finalize () called.\n");
      fflush(stderr);
    }
#endif // STEER
  
  if (net.id == 0)
    {
      total_time = myClock () - total_time;
      fprintf (timings_ptr, "total time (s):                            %.3f\n\n", total_time);

      fprintf (timings_ptr, "Sub-domains info:\n\n");
      
      for (int n = 0; n < net.procs; n++)
	{
	  fprintf (timings_ptr, "rank: %i, fluid sites: %i\n", n, proc_fluid_sites[ n ]);
	}
      
      fclose (timings_ptr);
    }
  free(proc_fluid_sites);
  
#ifndef NOMPI
  net.err = MPI_Finalize ();
#endif
  
  return(0);
}


