#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>

#define networkXDR

#ifdef networkXDR
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif


#define GLUTCALLBACK


#define PORT 65250

#ifdef MAC_GLUT
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "eVizRLEUtil.h"


#define GLUTCALLBACK
#define SCREEN_SIZE_MAX    1024 * 1024

int pixels_x, pixels_y;

unsigned char *pixel_data;
unsigned char *compressed_data;
int compressed_frame_size;
int frame_number;

int sockfd;

unsigned int width;
unsigned int height;
unsigned int bpp;

int pixels_max;

u_int sizeToRecv = (2 + SCREEN_SIZE_MAX) * sizeof(unsigned int);

char *xdrReceiveBuffer;


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
  int ret;
  int bytesReceived = 0;
  int lengthToReceive;
  
  XDR xdr_network_stream;
  
  unsigned short shortWidth, shortHeight;
  
  unsigned int four_compressed_data;
  
  int m, n;
  
  
  xdrmem_create (&xdr_network_stream, xdrReceiveBuffer, sizeToRecv, XDR_DECODE);
  
  printf ("%i\n", xdr_getpos (&xdr_network_stream));
  
  lengthToReceive = 4;
  recv_all (sockfd, xdrReceiveBuffer, &lengthToReceive);
  
  bytesReceived += lengthToReceive;
  
  lengthToReceive = 4;
  recv_all (sockfd, xdrReceiveBuffer, &lengthToReceive);
  
  bytesReceived += lengthToReceive;
  
  xdr_int (&xdr_network_stream, &compressed_frame_size);
  xdr_int (&xdr_network_stream, &frame_number);
  
  printf("got compressed frame size = %i, frame #%i\n", compressed_frame_size, frame_number);
  
  xdr_setpos (&xdr_network_stream, 0);
  
  lengthToReceive = compressed_frame_size;//*4;
  recv_all (sockfd, xdrReceiveBuffer, &lengthToReceive);
  bytesReceived += lengthToReceive;
  
  printf("received %i Bytes\n", bytesReceived);
  
  //for (int i = 0; i < compressed_frame_size; i++)
  //  {
  //    xdr_u_char (&xdr_network_stream, &compressed_data[i]);
  //  }
  m = (compressed_frame_size >> 2) << 2;
  
  if (m < compressed_frame_size) m += 4;
  
  n = 0;
  
  for (int i = 0; i < (m >> 2); i++)
    {
      xdr_u_int (&xdr_network_stream, &four_compressed_data);
      
      compressed_data[n++] = four_compressed_data & (1U << 8U);
      
      if (n++ < compressed_frame_size)
  	{
  	  compressed_data[n] = (four_compressed_data >> 8U) & (1U << 8U);
  	}
      if (n++ < compressed_frame_size)
  	{
  	  compressed_data[n] = (four_compressed_data >> 16U) & (1U << 8U);
  	}
      if (n++ < compressed_frame_size)
  	{
  	  compressed_data[n] = (four_compressed_data >> 24U) & (1U << 8U);
  	}
    }
  
  xdr_destroy (&xdr_network_stream);
  
  printf ("got frame # %i, # pixels %i\n", frame_number, compressed_frame_size);
  fflush (NULL);
  
  printf ("got the entire frame\n");
  fflush (NULL);
  
  ((unsigned char*)&shortWidth)[0] = compressed_data[0];
  ((unsigned char*)&shortWidth)[1] = compressed_data[1];
  ((unsigned char*)&shortHeight)[0] = compressed_data[2];
  ((unsigned char*)&shortHeight)[1] = compressed_data[3];
  
  width = shortWidth;
  height = shortHeight;
  
  pixels_x = width;
  pixels_y = height;
  
  pixels_max = width * height;
  
  if (pixel_data == NULL)
    {
      pixel_data = (unsigned char *)malloc(sizeof(unsigned char) * 3 * pixels_max);
    }
  
  eViz_RLE_readMemory (compressed_data, compressed_frame_size, &width, &height, &bpp, pixel_data);
  
  printf("size %i, width %i, height %i, bpp %i\n", compressed_frame_size, width, height, bpp);
  fflush(NULL);
}

void DisplayFrame ()
{
  int k;
  
  
  glPointSize (1.F);
  glBegin (GL_POINTS);
  
  for (int i = 0; i < width; i++)
    {
      for (int j = 0; j < height; j++)
	{
	  k = (i * height + j) * bpp;
	  
	  glColor3f (pixel_data[ k   ] * (1.F / 255.F),
		     pixel_data[ k+1 ] * (1.F / 255.F),
		     pixel_data[ k+2 ] * (1.F / 255.F));
	  
	  glVertex2f (i, j);
	}
    }
  glEnd();
  
  ReceiveFrame ();
}


void GLUTCALLBACK Display (void)
{
  glClear (GL_COLOR_BUFFER_BIT);
  
  DisplayFrame ();
  
  glutSwapBuffers ();
}


void GLUTCALLBACK Reshape (GLsizei w, GLsizei h)
{
  pixels_x = w;
  pixels_y = h;
  
  glViewport(0, 0, w, h);
  
  glLoadIdentity ();
  gluOrtho2D (0.F, (float)pixels_x, 0.F, (float)pixels_y);
  
  if (pixels_x * pixels_y > pixels_max)
    {
      pixels_max = pixels_x * pixels_y;
      pixel_data = (unsigned char *)realloc(pixel_data,
					    sizeof(unsigned char) * 3 * pixels_max);
    }
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


int main(int argc, char *argv[])
{  
  struct hostent *he;
  struct sockaddr_in their_addr;
  
  int compressed_frame_size_max;
  
  
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
  
  compressed_frame_size_max = 1024 * 1024;
  compressed_data = (unsigned char *)malloc(sizeof(unsigned char) * compressed_frame_size_max);
  
  pixel_data = NULL;
  
  xdrReceiveBuffer = (char *)malloc(sizeof(char) * sizeToRecv);
  
  ReceiveFrame ();

  glutInit (&argc, argv);
  OpenWindow (pixels_x, pixels_y);
  
  glLoadIdentity ();
  gluOrtho2D (0.F, (float)pixels_x, 0.F, (float)pixels_y);
  glClearColor (1.F, 1.F, 1.F, 0.F);
  
  glutReshapeFunc (Reshape);
  glutIdleFunc (Display);
  glutDisplayFunc (Display);
  glutMainLoop ();
  
  close (sockfd);
  
  free(xdrReceiveBuffer);
  free(compressed_data);
  free(pixel_data);
  
  return 0;
} 

