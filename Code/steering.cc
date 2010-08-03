#include "xdrWriter.h"
#include "xdrMemWriter.h"

#ifndef NO_STEER

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#include <netdb.h>

#include <sys/time.h>
#include <sched.h>

#include <sys/stat.h>
#include <string.h>

#include "network.h"
#include "steering.h"
#include "colourpalette.h"
#include "visthread.h"
#include "steering-sim-params.h"

#ifdef _AIX
#include <fcntl.h>
#endif

#include "http_post.h"

#define MYPORT 65250
#define CONNECTION_BACKLOG 10

pthread_mutex_t steer_param_lock = PTHREAD_MUTEX_INITIALIZER;

extern bool updated_mouse_coords;

char host_name[255];

float steer_par[ STEERABLE_PARAMETERS + 1 ] = {0.0, 0.0, 0.0,    // scene center (dx,dy,dz)
					       45.0, 45.0,       // longitude and latitude
					       1.0, 0.1,         // zoom and brightness
					       0.1, 0.1,         // velocity and stress ranges
					       80.0, 120.0,      // Minimum pressure and maximum pressure for Colour mapping
					       1.0,              // Glyph length
					       512, 512,         // Rendered frame size, pixel x and pixel y
					       -1.0, -1.0,         // x-y position of the mouse of the client
					       0.0,              // signal useful to terminate the simulation
					       0.0, 	         // Vis_mode 
					       5.0,	         // vis_streaklines_per_pulsatile_period
					       100.0,	         // vis_streakline_length					   
					       0.0};             // doRendering


double frameTiming()
{
  struct timeval time_data;
  gettimeofday (&time_data, NULL);
  return (double)time_data.tv_sec + (double)time_data.tv_usec / 1.0e6;
}

void *hemeLB_network (void *ptr)
{
  int steering_session_id = *(int*)ptr;
  char steering_session_id_char[255];
  
  
  snprintf(steering_session_id_char, 255, "%i", steering_session_id);
  
  setRenderState(0);
  
  gethostname (host_name, 255);
  
  FILE *f = fopen ("env_details.asc","w");
  
  fprintf (f, "%s\n", host_name);
  fclose (f);
  
  // fprintf (timings_ptr, "MPI 0 Hostname -> %s\n\n", host_name);
  
  //printf("kicking off network thread.....\n"); fflush(0x0);
  
  int sock_fd;
  int new_fd;
  int yes = 1;
  
  int is_broken_pipe = 0;
  int frame_number = 0;
  
  pthread_t steering_thread;
  pthread_attr_t steering_thread_attrib; 
  
  static char ip_addr[16];
  static char rank_0_host_details[1024];
  
  
  pthread_attr_init (&steering_thread_attrib);
  pthread_attr_setdetachstate (&steering_thread_attrib, PTHREAD_CREATE_JOINABLE);
  
  signal(SIGPIPE, SIG_IGN); // Ignore a broken pipe 
  
  HTTP::get_host_details(rank_0_host_details, ip_addr);
  
  HTTP::request("bunsen.chem.ucl.ac.uk", 28080, "/ahe/test/rendezvous/", steering_session_id_char, rank_0_host_details);
  
  while (1)
    {
      setRenderState(0);
      
      // pthread_mutex_lock (&LOCK);
      // sem_wait( &nrl );
      
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
      
      // socklen_t _length = sizeof(my_address);
      // 
      // getsockname (sock_fd, (struct sockaddr *) &my_address,&_length);
      // printf("Server Port is: %d\n", ntohs(my_address.sin_port));
      
      if ((new_fd = accept (sock_fd, (struct sockaddr *)&their_addr, &sin_size)) == -1)
	{
	  perror("accept");
	  // continue;
	}
      
      // fprintf (timings_ptr, "server: got connection from %s (FD %i)\n", inet_ntoa (their_addr.sin_addr), new_fd);
      // printf ("RG thread: server: got connection from %s (FD %i)\n", inet_ntoa (their_addr.sin_addr), new_fd);
      
      pthread_create (&steering_thread, &steering_thread_attrib, hemeLB_steer, (void*)&new_fd);	  
      
      close(sock_fd);
      
      is_broken_pipe = 0;
      
      // pthread_mutex_unlock ( &LOCK );
      sem_wait(&connected_sem); // grab the semaphore
      connected = 1;
      sem_post(&connected_sem);
      
      // setRenderState(1);
      // At this point we're ready to send a frame...
      // setRendering=1;
      // sem_wait( &nrl );
      
      while (!is_broken_pipe)
	{
	  // printf("THREAD: waiting for signal that frame is ready to send..\n"); fflush(0x0);
	  // pthread_mutex_lock ( &LOCK );
	  
	  bool is_frame_ready_local = 0;
	  
	  while(!is_frame_ready_local) {
	    usleep(5000);
	    sem_wait( &nrl );
	    is_frame_ready_local = is_frame_ready;
	    sem_post( &nrl );
	    // printf("THREAD is_frame_ready_local %i\n", is_frame_ready_local);
	  }
	  sem_wait( &nrl );
	  sending_frame = 1;
	  // printf("THREAD sending frame = 1\n");
	  // pthread_cond_wait (&network_send_frame, &LOCK);
	  // setRenderState(0);
	  // printf("THREAD: received signal that frame is ready to send..\n"); fflush(0x0);
	  
	  double frameTimeStart = frameTiming();		
	  
	  int bytesSent = 0;
	  
	int pixeldatabytes = 8;
	char xdr_pixel[pixeldatabytes];
	XdrMemWriter pixelWriter = XdrMemWriter(xdr_pixel, pixeldatabytes);
	
	pixelWriter.write(screen.pixels_x);
	pixelWriter.write(screen.pixels_y);

        
	Network::send_all(new_fd, xdr_pixel, &pixeldatabytes);

        XdrMemWriter pixelDataWriter = XdrMemWriter(xdrSendBuffer_pixel_data, pixel_data_bytes);
        XdrMemWriter frameDetailsWriter = XdrMemWriter(xdrSendBuffer_frame_details, frame_details_bytes);
	
	for (int i = 0; i < col_pixels_recv[RECV_BUFFER_A]; i++)
	{
          pixelDataWriter.writePixel (&col_pixel_recv[RECV_BUFFER_A][i], ColourPalette::PickColour);
	}
	
	int frameBytes = pixelDataWriter.getCurrentStreamPosition();
	
	frameDetailsWriter.write(frameBytes);
	
	int detailsBytes = frameDetailsWriter.getCurrentStreamPosition();
	
	int ret = Network::send_all(new_fd, xdrSendBuffer_frame_details, &detailsBytes);
	
	if (ret < 0) {
	  printf("RG thread: broken network pipe...\n");
	  is_broken_pipe = 1;
	  // pthread_mutex_unlock ( &LOCK );
	  sem_post(&nrl);
	  setRenderState(0);
	  break;
	} else {
	  bytesSent += detailsBytes;
	}
	
	ret = Network::send_all(new_fd, xdrSendBuffer_pixel_data, &frameBytes);
	
	if (ret < 0) {
	  printf("RG thread: broken network pipe...\n");
	  is_broken_pipe = 1;
	  // pthread_mutex_unlock ( &LOCK );
	  sem_post(&nrl);
	  setRenderState(0);
	  break;
	} else {
	  bytesSent += frameBytes;
	}
	
	simulationParameters* Sim = new simulationParameters();
	Sim->collectGlobalVals();
	int sizeToSend = Sim->sim_params_bytes;
	Network::send_all(new_fd, Sim->pack(), &sizeToSend);
	// printf ("Sim bytes sent %i\n", sizeToSend);
	delete Sim;
	
	// fprintf (timings_ptr, "bytes sent %i\n", bytesSent);
	// printf ("RG thread: bytes sent %i\n", bytesSent);
	
	setRenderState(1);
	
        double frameTimeSend = frameTiming() - frameTimeStart;
	
	// printf("Time to send frame = %0.6f s\n", frameTimeSend);
	
        double timeDiff = (1.0/25.0) - frameTimeSend;
	
        if (timeDiff > 0.0)
	  {
	    // printf("Sleeping for %0.6f s\n", timeDiff);
	    
	    usleep(timeDiff*1.0e6);
	  }
	  
	// pthread_mutex_unlock ( &LOCK );
	// sem_post(&nrl);
	
	sending_frame = 0;
	is_frame_ready = 0;
	sem_post( &nrl );
	
	frame_number++;
	
	} // while (is_broken_pipe == 0)
      
      close(new_fd);
      
      sem_wait(&connected_sem); // grab the semaphore
      connected = 0;
      sem_post(&connected_sem);
      
      // pthread_join(steering_thread, NULL);
      
    } // while(1)
}

void* hemeLB_steer (void* ptr)
{
  int read_fd = *(int*)ptr;
  
  // printf("Kicking off steering thread with FD %i\n", read_fd);
  
  int num_chars = STEERABLE_PARAMETERS * sizeof(float) / sizeof(char);
  int bytes = sizeof(char) * num_chars;
  
  char* xdr_steering_data = (char*)malloc(bytes);
  
  while(1)
    {
      XDR xdr_steering_stream;
    
      xdrmem_create(&xdr_steering_stream, xdr_steering_data, bytes, XDR_DECODE);
      
      while(1)
	{
	  struct timeval tv;
	  fd_set readfds;
	  
	  tv.tv_sec = 0;
	  tv.tv_usec = 0;
	  
	  FD_ZERO(&readfds);
	  FD_SET(read_fd, &readfds);
	  
	  select(read_fd+1, &readfds, NULL, NULL, &tv);
	  // printf("STEERING: Polling..\n"); fflush(0x0);  
	  
	  int ret;	
	  
	  if(FD_ISSET(read_fd, &readfds))
	    {
	      /* If there's something to read, read it... */
	      //				printf("STEERING: Got data\n"); fflush(0x0);
	      ret = Network::recv_all(read_fd, xdr_steering_data, &num_chars);
	      sched_yield();
	      break;
	    } else {
	    usleep(5000);	
	  }
	  
	  if (ret < 0)
	    {
	      printf("Steering thread: broken network pipe...\n");
	      free(xdr_steering_data);
	      xdr_destroy(&xdr_steering_stream);
	      return NULL;
	    }
	  sched_yield();
	}
      // pthread_mutex_lock(&steer_param_lock);
      
      sem_wait(&steering_var_lock);
      
      for (int i = 0; i < STEERABLE_PARAMETERS; i++)
	xdr_float(&xdr_steering_stream, &steer_par[i]);
      
      sem_post(&steering_var_lock);
      
      /* printf("Got steering params ");
	 for (int i = 0; i < STEERABLE_PARAMETERS; i++) 
	 printf("%0.4f ", steer_par[i]);
	 printf("\n"); */
      
      if (steer_par[14] > -1.0 && steer_par[15] > -1.0)
	updated_mouse_coords = 1;  
      
      // pthread_mutex_unlock(&steer_param_lock);
      
      xdr_destroy(&xdr_steering_stream);
    }
  free(xdr_steering_data);
  
  return NULL;
}


#else // NO_STEER

float steer_par[ STEERABLE_PARAMETERS + 1 ] = {0.0, 0.0, 0.0,    // scene center (dx,dy,dz)
					       45.0, 45.0,          // longitude and latitude
					       1.0, 0.03,        // zoom and brightness
					       0.1, 0.1,         // velocity and stress ranges
					       80.0, 120.0,      // Minimum pressure and maximum pressure for Colour mapping
					       1.0,              // Glyph length
					       512, 512,         // Rendered frame size, pixel x and pixel y
					       -1.0, -1.0,         // x-y position of the mouse of the client
					       0.0,              // signal useful to terminate the simulation
					       0.0, 	         // Vis_mode 
					       5.0,	         // vis_streaklines_per_pulsatile_period
					       100.0,	         // vis_streakline_length					   
					       0.0};             // doRendering
#endif // NO_STEER
