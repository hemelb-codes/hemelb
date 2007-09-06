#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>

#define GLUTCALLBACK

#define PORT 10000

#include <GL/glut.h>

#include "eVizRLEUtil.h"


int pixels_x = 512, pixels_y = 512;

int frame_number;
unsigned char* pixel_data;
unsigned char* compressed_data;
int compressed_frame_size;

int sockfd;

unsigned int width;
unsigned int height;
unsigned int bpp;


void ReceiveAndFillDispBuffers ()
{
  int pixel_i, pixel_j;
  int i, k;
  
  
  recv(sockfd, &frame_number, sizeof(frame_number), 0);
  recv(sockfd, &compressed_frame_size, sizeof(compressed_frame_size), 0);
  
  printf("got frame # %i, # pixels %i\n", frame_number, compressed_frame_size);
  fflush(NULL);
  
  for (i = 0; i < compressed_frame_size; i++)
    {
      recv(sockfd, &compressed_data[i], sizeof(compressed_data[i]), 0);
    } 
  
  printf("got the entire frame\n");
  fflush(NULL);
  
  eViz_RLE_readMemory (compressed_data, compressed_frame_size, &width, &height, &bpp, pixel_data);
  
  printf("size %i, width %i, height %i, bpp %i\n", compressed_frame_size, width, height, bpp);
  
  glPointSize (1.F);
  
  glBegin (GL_POINTS);
  
  pixel_i = 0;
  pixel_j = 0;
  
  for (i = 0; i < width * height; i++)
    {
      if (i > 0 && i%width == 0)
	{
	  pixel_i = 0;
	  pixel_j++;
	}
      k = (pixel_i * width + pixel_j) * bpp;
      
//      float r, g, b;
//      memcpy(&r, &pixel_data[ k   ], 4);
//      memcpy(&g, &pixel_data[ k+4 ], 4);
//      memcpy(&b, &pixel_data[ k+8 ], 4);
//      
//      if (r == 0.F && g == 0.F && b == 0.F)
//	{
//	  r = 1.F;
//	  g = 1.F;
//	  b = 1.F;
//	}

      glColor3f(pixel_data[ k   ] * (1.F / 255.F),
		pixel_data[ k+1 ] * (1.F / 255.F),
		pixel_data[ k+2 ] * (1.F / 255.F));
      
      glVertex2f (pixel_i, pixel_j);
      
      ++pixel_i;
    }
  glEnd();
}


void GLUTCALLBACK Display (void)
{
  glClear (GL_COLOR_BUFFER_BIT);
  ReceiveAndFillDispBuffers ();
  glutSwapBuffers ();
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
  pixel_data = (unsigned char *)malloc(sizeof(unsigned char) * 3 * pixels_x * pixels_y);
  compressed_data = (unsigned char *)malloc(sizeof(unsigned char) * 3 * pixels_x * pixels_y);
  
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
  
  if (connect (sockfd, (struct sockaddr *)&their_addr, sizeof their_addr) == -1)
    {
      perror("connect");
      exit(1);
    }
  
  glutInit (&argc, argv);
  OpenWindow (pixels_x, pixels_y);
  
  glLoadIdentity ();
  gluOrtho2D (0.F, 512.F, 0.F, 512.F);
  glClearColor (1.F, 1.F, 1.F, 0.F);
  
  glutIdleFunc (Display);
  glutDisplayFunc (Display);
  glutMainLoop ();
  
  close (sockfd);
  
  return 0;
} 

