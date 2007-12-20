#ifdef RG

#include "config.h"
#include "network.h"

#include <pthread.h>

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

char host_name[255];

pthread_mutex_t network_buffer_copy_lock;
pthread_cond_t network_send_frame;
int send_array_length;

pthread_t network_thread;
pthread_attr_t pthread_attrib;

char *xdrSendBuffer_pixel_data;
char *xdrSendBuffer_frame_details;

// data per pixel are colour id and pixel id (2 * sizeof(int) bytes)
int data_per_pixel = 2;
int bytes_per_pixel_data = data_per_pixel * sizeof(int);

// one int for colour_id and one for pixel id
u_int pixel_data_bytes = IMAGE_SIZE * bytes_per_pixel_data;

// it is assumed that the frame size is the only detail
u_int frame_details_bytes = 1 * sizeof(int);

int bits_per_char = sizeof(char) * 8;
int bits_per_two_chars = 2 * bits_per_char;

void setupNetworkBuffersAndThread() {
	
	xdrSendBuffer_pixel_data = (char *)malloc(pixel_data_bytes);
	xdrSendBuffer_frame_details = (char *)malloc(frame_details_bytes);
	
	pthread_mutex_init (&network_buffer_copy_lock, NULL);
	pthread_cond_init (&network_send_frame, NULL);
	
	pthread_attr_init (&pthread_attrib);
	pthread_attr_setdetachstate (&pthread_attrib, PTHREAD_CREATE_JOINABLE);
	
	pthread_create (&network_thread, &pthread_attrib, hemeLB_network, NULL);
	
}

void CleanUpNetworkBuffersAndThread() {
		// there are some problems if the following function is called
		//pthread_join (network_thread, NULL);
		free(xdrSendBuffer_frame_details);
		free(xdrSendBuffer_pixel_data);
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

void *hemeLB_network (void *ptr) {
	
	gethostname (host_name, 255);
	
	FILE *f = fopen ("env_details.asc","w");
	
	fprintf (f, "%s\n", host_name);
	fclose (f);
	
	// fprintf (timings_ptr, "MPI 0 Hostname -> %s\n\n", host_name);
	
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
	
	while (1) {
		
		pthread_mutex_lock ( &network_buffer_copy_lock );
		
		struct sockaddr_in my_address;
		struct sockaddr_in their_addr; // client address
		
		socklen_t sin_size;
		
		
		if ((sock_fd = socket (AF_INET, SOCK_STREAM, 0)) == -1) {
			perror("socket");
			exit (1);
		}
		
		if (setsockopt (sock_fd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1) {
			perror("setsockopt");
			exit (1);
		}
		
		my_address.sin_family = AF_INET;
		my_address.sin_port = htons (MYPORT);
		my_address.sin_addr.s_addr = INADDR_ANY;
		memset (my_address.sin_zero, '\0', sizeof my_address.sin_zero);
		
		if (bind (sock_fd, (struct sockaddr *)&my_address, sizeof my_address) == -1) {
			perror ("bind");
			exit (1);
		}
		
		if (listen (sock_fd, CONNECTION_BACKLOG) == -1) {
			perror ("listen");
			exit (1);
		}
		
		sin_size = sizeof their_addr;
		
		if ((new_fd = accept (sock_fd, (struct sockaddr *)&their_addr, &sin_size)) == -1) {
			perror("accept");
			continue;
		}
		
//		fprintf (timings_ptr, "server: got connection from %s\n", inet_ntoa (their_addr.sin_addr));
		printf ("server: got connection from %s\n", inet_ntoa (their_addr.sin_addr));
		
		close(sock_fd);
		
		is_broken_pipe = 0;
		
		pthread_mutex_unlock ( &network_buffer_copy_lock );
		
		while (!is_broken_pipe) {
			pthread_mutex_lock ( &network_buffer_copy_lock );
			pthread_cond_wait (&network_send_frame, &network_buffer_copy_lock);
			
			int bytesSent = 0;
			
			XDR xdr_network_stream_frame_details;
			XDR xdr_network_stream_pixel_data;
			
			xdrmem_create (&xdr_network_stream_pixel_data, xdrSendBuffer_pixel_data,
						   pixel_data_bytes, XDR_ENCODE);
			
			xdrmem_create (&xdr_network_stream_frame_details, xdrSendBuffer_frame_details,
						   frame_details_bytes, XDR_ENCODE);
			
			for (int i = 0; i < vis.col_pixels_locked; i++) {
				
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
			
//			fprintf (timings_ptr, "bytes sent %i\n", bytesSent);
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