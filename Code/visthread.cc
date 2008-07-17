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

#include <sys/stat.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include <string.h>

#include "network.h"
#include "config.h"
#include "colourpalette.h"
#include "steering.h"
#include "visthread.h"

#define MYPORT 65250
#define CONNECTION_BACKLOG 10

int doRendering = 0;
int ShouldIRenderNow = 0;

pthread_mutex_t var_lock = PTHREAD_MUTEX_INITIALIZER;

// data per pixel inlude the data for the pixel location and 4 colours
//in RGB format, thus #bytes per pixel are (sizeof(int)+4
//rgb)=sizeof(int)+4*3*sizeof(unsigned char))
int bytes_per_pixel_data = sizeof(int) + 4 * sizeof(unsigned char);

// one int for colour_id and one for pixel id
u_int pixel_data_bytes = IMAGE_SIZE * bytes_per_pixel_data;

// it is assumed that the frame size is the only detail
u_int frame_details_bytes = 1 * sizeof(int);
char* xdrSendBuffer_pixel_data;
char* xdrSendBuffer_frame_details;

char host_name[255];

void setRenderState(int val) {
  pthread_mutex_lock(&var_lock);
//  doRendering = val;
  ShouldIRenderNow = val;
  pthread_mutex_unlock(&var_lock);
}

void *hemeLB_network (void *ptr) {

  setRenderState(0);

  gethostname (host_name, 255);
  
  FILE *f = fopen ("env_details.asc","w");
  
  fprintf (f, "%s\n", host_name);
  fclose (f);
  
  // fprintf (timings_ptr, "MPI 0 Hostname -> %s\n\n", host_name);

  printf("kicking off network thread.....\n"); fflush(0x0);
  
  int sock_fd;
  int new_fd;
  int yes = 1;
  
  int is_broken_pipe = 0;
  int frame_number = 0;
  
  pthread_t steering_thread;
  pthread_attr_t steering_thread_attrib; 
  pthread_attr_init (&steering_thread_attrib);
  pthread_attr_setdetachstate (&steering_thread_attrib, PTHREAD_CREATE_JOINABLE);

  signal(SIGPIPE, SIG_IGN); // Ignore a broken pipe 
  
  while (1)
    {

      setRenderState(0);

      pthread_mutex_lock (&LOCK);
	    
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
      
      pthread_mutex_unlock ( &LOCK );

      setRenderState(1);

	// At this point we're ready to send a frame...
	
//	setRendering=1;
      
     while (!is_broken_pipe)
	 {

	  printf("THREAD: waiting for signal that frame is ready to send..\n"); fflush(0x0);

	  pthread_mutex_lock ( &LOCK );
	  pthread_cond_wait (&network_send_frame, &LOCK);

      setRenderState(0);

	  printf("THREAD: received signal that frame is ready to send..\n"); fflush(0x0);

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
		    printf("RG thread: broken network pipe...\n");
            is_broken_pipe = 1;
			pthread_mutex_unlock ( &LOCK );
            setRenderState(0);
            break;
          } else {
            bytesSent += detailsBytes;
          }
	  
	  ret = send_all(new_fd, xdrSendBuffer_pixel_data, &frameBytes);
	  
          if (ret < 0) {
		    printf("RG thread: broken network pipe...\n");
            is_broken_pipe = 1;
			pthread_mutex_unlock ( &LOCK );
            setRenderState(0);
            break;
          } else {
            bytesSent += frameBytes;
          }
	  
	  //fprintf (timings_ptr, "bytes sent %i\n", bytesSent);
	  printf ("RG thread: bytes sent %i\n", bytesSent);

      setRenderState(1);
	  
      xdr_destroy (&xdr_network_stream_frame_details);
	  xdr_destroy (&xdr_network_stream_pixel_data);

	  pthread_mutex_unlock ( &LOCK );
	  
	  frame_number++;
	  
	} // while (is_broken_pipe == 0)
      
      close(new_fd);
      
    } // while(1)
}

