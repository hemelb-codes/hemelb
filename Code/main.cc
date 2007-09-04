// In this file, the functions useful to initiate/end the LB simulation
// and perform the dynamics are reported

#include "config.h"

RT rt;

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

#define MYPORT 10000

Pixel* send_array;
int send_array_length;
int send_frame_count = 0;

pthread_mutex_t network_buffer_copy_lock;
pthread_cond_t network_send_frame;

#endif // RG


#ifdef RG

void *hemeLB_network(void *ptr) {

	printf("kicking off thread.....\n");

	while(true) {

	int sock_fd, new_fd;
	int yes = 1;
	struct sockaddr_in my_address;
	struct sockaddr_in their_addr; // client address

	socklen_t sin_size;

	if ((sock_fd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
	perror("socket");
	exit(1);
	}

	if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1) {
		perror("setsockopt");
		exit(1);
	}

	my_address.sin_family = AF_INET;
	my_address.sin_port = htons(MYPORT);
	my_address.sin_addr.s_addr = INADDR_ANY;
	memset(my_address.sin_zero, '\0', sizeof my_address.sin_zero);

	if (bind(sock_fd, (struct sockaddr *)&my_address, sizeof my_address) == -1) {
		perror("bind");
		exit(1);
	}

	if (listen(sock_fd, 10) == -1) {
		perror("listen");
		exit(1);
	}

	sin_size = sizeof their_addr;

	if ((new_fd = accept(sock_fd, (struct sockaddr *)&their_addr, &sin_size)) == -1)
		perror("accept");

	printf("server: got connection from %s\n", inet_ntoa(their_addr.sin_addr));

	close(sock_fd);

	int frameNumber = 0;

	while(true) {

		pthread_mutex_lock( &network_buffer_copy_lock );
		pthread_cond_wait(&network_send_frame, &network_buffer_copy_lock);

		send(new_fd, &frameNumber, sizeof(frameNumber), 0);

		send(new_fd, &compressedFrameSize, sizeof(compressedFrameSize), 0);

		for(int i=0; i<compressedFrameSize; i++)
			send(new_fd, &compressedData[i], sizeof(compressedData[i]), 0);

		printf("done sending array...");

		pthread_mutex_unlock( &network_buffer_copy_lock );

		frameNumber++;

	}

	close(new_fd);

	}

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
  printf ("Usage: %s input file name with config, parameters, output,\n", progname);
  printf ("          checkpoint, rt_parameters and image file names\n");
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
  
  char *input_file_name(argv[1]);
  char input_config_name[800];
  char input_parameters_name[800];
  char output_config_name[800];
  char checkpoint_config_name[800];
  char rt_parameters_name[800];
  char output_image_name[800];
  
  FILE *input_file = fopen (input_file_name, "r");

  fscanf (input_file, "%s ", input_config_name);
  fscanf (input_file, "%s ", input_parameters_name);
  fscanf (input_file, "%s ", output_config_name);
  fscanf (input_file, "%s ", checkpoint_config_name);
  fscanf (input_file, "%s ", rt_parameters_name);
  fscanf (input_file, "%s ", output_image_name);
  
  fclose (input_file);
  
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
  
  lbmReadParameters (input_parameters_name, &lbm, &net);

#ifdef RG

	if(net.id == 0) {

		pthread_mutex_init(&network_buffer_copy_lock, NULL);
		pthread_cond_init(&network_send_frame, NULL);

		pthread_attr_init(&pthread_attrib);
		pthread_attr_setdetachstate(&pthread_attrib, PTHREAD_CREATE_JOINABLE);

		pthread_create(&network_thread, &pthread_attrib, hemeLB_network, NULL);

	}

#endif // RG

  lbmInit (input_config_name, checkpoint_config_name, &lbm, &net);
  
  netInit (&lbm, &net, &rt);
  
  rtReadParameters (rt_parameters_name, &rt, &net);
  
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
  
  for (time_step = 1; time_step <= lbm.time_steps_max; time_step++)
    {
      write_checkpoint = 0;
      check_convergence = 0;
      perform_rt = 0;

	if(net.id==0) { printf("time step %i\n",time_step); fflush(0x0); }
      
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
      stability = lbmCycle (write_checkpoint, check_convergence, perform_rt,
			    &is_converged, &lbm, &net);
      
      if (perform_rt)
	{
	  rtRayTracing (AbsorptionCoefficients, &lbm, &net, &rt);
	}
      
      if (stability == UNSTABLE || is_converged) break;
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
      printf ("1) rank\n");
      printf ("2) collisions at interface-dependent lattice sites\n");
      printf ("3) copying to the buffers to send/recv\n");
      printf ("4) communicational time\n");
      printf ("5) collisions + streaming from inner sites\n");
      printf ("6) streaming from interface-dependent lattice sites\n");
      printf ("7) all-reduce + convergence test + checkpoint\n");
      printf ("8) ray tracing \n");
      printf ("9) output (1 time only, seconds) \n");
      fflush (stdout);
    }
  total_time = 0.;
  
  for (n = 0; n < 7; n++) total_time += net.timing[ n ];
  
  for (n = 0; n < 7; n++)
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
	  else if (n == 6)
	    {
	      net.timing[ n ] /= rt.ray_tracing_count;
	    }
	}
    }
  
  printf ("%i %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e 1) -> 9)\n", net.id,
	  net.timing[0], net.timing[1], net.timing[2], net.timing[3],
	  net.timing[4], net.timing[5], net.timing[6], net.timing[7]);
  fflush (stdout);
  
  rtEnd (&rt);
  lbmEnd (&lbm);
  netEnd (&net, &rt);

#ifdef RG

	pthread_join(network_thread, NULL);

#endif // RG

  net.err = MPI_Finalize ();
  
  return(0);
}

