#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <netdb.h>
#include <time.h>
#include <math.h>

//#ifndef XT3
//#include <sys/types.h>
//#include <netinet/in.h>
//#include <sys/socket.h>
//#endif

#ifdef XT3
#include "types.h"
#include "xdr.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#ifndef GLUTCALLBACK
#define GLUTCALLBACK
#endif

#define PORT 65250

#ifdef MAC_GLUT
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#define PIXELS_X      1024
#define PIXELS_Y      1024
#define IMAGE_SIZE    PIXELS_X * PIXELS_Y


struct PixelData
{
  unsigned char r, g, b;
  
  short int i, j;
};

PixelData pixel_data[ IMAGE_SIZE ];


int pixels_x, pixels_y;

int col_pixels;
int frame_size;
int frame_number;

int sockfd;

unsigned int width;
unsigned int height;
unsigned int bpp;

// data per pixel are colour id and pixel id (2 * sizeof(int) * bytes)
int data_per_pixel = 2;
int bytes_per_pixel_data = data_per_pixel * sizeof(int);

int bits_per_char = sizeof(char) * 8;
int bits_per_two_chars = 2 * bits_per_char;

// r-, g-, b- pixel components need "bits_per_char" bits
int colour_mask = (1 << bits_per_char) - 1;

// x-, y- pixel components need "bits_per_two_chars" bits
int pixel_mask = (1 << bits_per_two_chars) - 1;

u_int pixel_data_bytes = IMAGE_SIZE * bytes_per_pixel_data;

// it is assumed that the frame size is the only detail
u_int frame_details_bytes = 1 * sizeof(int);

char *xdrReceiveBuffer;
char *xdrReceiveBuffer_size;


double myClock ()
{
  struct timeval time_data;
  gettimeofday (&time_data, NULL);
  return (double)time_data.tv_sec + (double)time_data.tv_usec / 1.e6;
}


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


void ReceiveFrame ()
{
  double start_time, transfer_time;
  
  int bytesReceived = 0;
  int detailsBytes, frameBytes;
  int colour_id, pixel_id;
  
  PixelData *pixel_data_p;
  
  
  XDR xdr_network_stream_size;
  XDR xdr_network_stream;
  
  
  xdrmem_create (&xdr_network_stream, xdrReceiveBuffer, pixel_data_bytes, XDR_DECODE);
  xdrmem_create (&xdr_network_stream_size, xdrReceiveBuffer_size, frame_details_bytes, XDR_DECODE);
  
  // printf ("getpos: %i\n", xdr_getpos (&xdr_network_stream));
  
  start_time = myClock();
  
  detailsBytes = frame_details_bytes;
  recv_all (sockfd, xdrReceiveBuffer_size, &detailsBytes);
  bytesReceived += detailsBytes;
  
  xdr_int (&xdr_network_stream_size, &frame_size);
  
  printf("got frame size = %i bytes\n", frame_size);
  // xdr_setpos (&xdr_network_stream, 0);
  
  frameBytes = frame_size;
  
  recv_all (sockfd, xdrReceiveBuffer, &frameBytes);
  
  transfer_time = myClock() - start_time;
  
  bytesReceived += frameBytes;
  
  printf("received %i Bytes in %0.4fs (%0.2fKB/s)\n", bytesReceived, transfer_time, bytesReceived/(transfer_time*1024.));
  
  col_pixels = frame_size / bytes_per_pixel_data;
  
  for (int i = 0; i < col_pixels; i++)
    {
      xdr_int (&xdr_network_stream, &colour_id);
      xdr_int (&xdr_network_stream, &pixel_id);
      
      pixel_data_p = &pixel_data[ i ];
      
      pixel_data_p->r = (colour_id >> bits_per_two_chars) & colour_mask;
      pixel_data_p->g = (colour_id >> bits_per_char     ) & colour_mask;
      pixel_data_p->b = (colour_id                      ) & colour_mask;
      
      pixel_data_p->i = (pixel_id >> bits_per_two_chars) & pixel_mask;
      pixel_data_p->j = (pixel_id                      ) & pixel_mask;
    }
  
  xdr_destroy (&xdr_network_stream);
  xdr_destroy (&xdr_network_stream_size);
  
  printf ("got the entire frame\n");
  fflush (NULL);
}


void DisplayFrame ()
{
  float scale_x = 1.F / (float)PIXELS_X;
  float scale_y = 1.F / (float)PIXELS_Y;
  
  PixelData *pixel_data_p;
  
  
  ReceiveFrame ();
  
  glPointSize (1.F);
  glBegin (GL_POINTS);
  
  for (int i = 0; i < col_pixels; i++)
    {
      pixel_data_p = &pixel_data[ i ];
      
      glColor3f (pixel_data_p->r * (1.F / 255.F),
		 pixel_data_p->g * (1.F / 255.F),
		 pixel_data_p->b * (1.F / 255.F));
      
      glVertex2f (-0.5F + scale_x * pixel_data_p->i,
		  -0.5F + scale_y * pixel_data_p->j);
    }
  glEnd ();
  glFlush ();
}


void SaveWindowImage (int pixels_x, int pixels_y, char *file_name)
{
  FILE *ppm_image_file_ptr = fopen (file_name, "wb");
  
  int i, j;
  
  unsigned char *data = NULL;  
  unsigned char *data_p = NULL;
  unsigned char *buffer;
  
  
  glReadBuffer (GL_FRONT);
  
  data = (unsigned char *)malloc(sizeof(unsigned char) * pixels_x * pixels_y * 3);
  
  buffer = (unsigned char *)malloc(sizeof(unsigned char) * pixels_x * 3);
  
  data_p = data;

  for (j = 0; j < pixels_y; j++)
    {
      glReadPixels (0, j, pixels_x, 1, GL_RGB, GL_UNSIGNED_BYTE, buffer);
      
      for (i = 0; i < pixels_x; i++)
	{
	  *data_p = buffer[ i * 3     ]; data_p++;
	  *data_p = buffer[ i * 3 + 1 ]; data_p++;
	  *data_p = buffer[ i * 3 + 2 ]; data_p++;
	}
    }
  
  free((unsigned char *)buffer);
  
  fprintf (ppm_image_file_ptr, "P6\n%i %i\n255\n", pixels_x, pixels_y);
  
  for (j = pixels_y - 1; j >= 0; j--)
    {
      fwrite (data + j * pixels_x * 3, 1, pixels_x * 3, ppm_image_file_ptr);
    }
  
  free((unsigned char *)data);
  
  fclose (ppm_image_file_ptr);
}


void GLUTCALLBACK Display (void)
{
  glClear (GL_COLOR_BUFFER_BIT);
  
  DisplayFrame ();
  
  glutSwapBuffers ();
}


void GLUTCALLBACK Reshape (GLsizei w, GLsizei h)
{
  //// the window is reshaped if necessary
  
  float ortho_x = 0.5F * (float)w / (float)PIXELS_X;
  float ortho_y = 0.5F * (float)h / (float)PIXELS_Y;
  
  pixels_x = w;
  pixels_y = h;
  
  glViewport(0, 0, w, h);
  
  glLoadIdentity ();
  gluOrtho2D (-ortho_x, ortho_x, -ortho_y, ortho_y);
}


void OpenWindow (int pixels_x, int pixels_y)
{
  glutInitDisplayMode (GLUT_RGB | GLUT_DOUBLE);
  glutInitWindowPosition (0, 0);
  glutInitWindowSize (pixels_x, pixels_y);
  
  glutCreateWindow (" ");
  
  glDisable (GL_DEPTH_TEST);
  glDisable (GL_BLEND);
  glShadeModel (GL_FLAT);
  glDisable (GL_DITHER);
  
  glClear (GL_COLOR_BUFFER_BIT);
}


void GLUTCALLBACK KeybordFunction (unsigned char key, int x, int y)
{
  if (key == 'S')
    {
      SaveWindowImage (pixels_x, pixels_y, "./image.ppm");
    }
  else if (key == 'q')
    {
      close (sockfd);
      
      free(xdrReceiveBuffer_size);
      free(xdrReceiveBuffer);
      
      exit(0);
    }
}


int main(int argc, char *argv[])
{  
  struct hostent *he;
  struct sockaddr_in their_addr;
  
  
  if (argc != 2)
    {
      fprintf(stderr, "usage: client hostname\n");
      exit(1);
    }
  
  if ((he = gethostbyname(argv[1])) == NULL)
    {
      herror("gethostbyname");
      exit(1);
    }
  
  if ((sockfd = socket(PF_INET, SOCK_STREAM, 0)) == -1)
    {
      perror("socket");
      exit(1);
    }
  
  their_addr.sin_family = AF_INET;
  their_addr.sin_port = htons(PORT);
  their_addr.sin_addr = *((struct in_addr *)he->h_addr);
  memset (their_addr.sin_zero, '\0', sizeof their_addr.sin_zero);
  
  if (connect (sockfd, (struct sockaddr *)&their_addr, sizeof their_addr) == -1) {
    perror("connect");
    exit(1);
  }
  
  
  pixels_x = 1024;
  pixels_y = 1024;
  
  xdrReceiveBuffer = (char *)malloc(pixel_data_bytes);
  xdrReceiveBuffer_size = (char *)malloc(frame_details_bytes);
  
  glutInit (&argc, argv);
  OpenWindow (pixels_x, pixels_y);
  
  glLoadIdentity ();
  gluOrtho2D (-0.5, 0.5F, -0.5F, 0.5F);
  
  glClearColor (1.F, 1.F, 1.F, 0.F);
  
  glutReshapeFunc (Reshape);
  glutIdleFunc (Display);
  glutDisplayFunc (Display);
  glutKeyboardFunc (KeybordFunction);
  glutMainLoop ();
  
  return 0;
}
