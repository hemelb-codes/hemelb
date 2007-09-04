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

int frameNumber;
unsigned char* pixelData;
unsigned char* compressedData;
int compressedFrameSize;

int sockfd;

void ReceiveAndFillDispBuffers() {

	recv(sockfd, &frameNumber, sizeof(frameNumber),0);
	recv(sockfd, &compressedFrameSize, sizeof(compressedFrameSize),0);

	printf("got frame # %i, # pixels %i\n", frameNumber, compressedFrameSize); fflush(NULL);

	for(int i=0; i<compressedFrameSize; i++) {
		recv(sockfd, &compressedData[i], sizeof(compressedData[i]), 0);
	} 

	printf("got the entire frame\n"); fflush(NULL);

	unsigned int width;
	unsigned int height;
	unsigned int bpp;

	eViz_RLE_readMemory(compressedData, compressedFrameSize, &width, &height, &bpp, pixelData);

	printf("size %i, width %i, height %i, bpp %i\n", compressedFrameSize, width, height, bpp);

	glPointSize (1.F);

	glBegin (GL_POINTS);

//	float r, g, b;

	int x_coord=0, y_coord=0;

	for(int i=0; i<width*height; i++) {

		if( i>0 && i%width==0 ) {
			x_coord = 0;
			y_coord++;
		}

		int I = (x_coord + y_coord*width)*bpp;

//		memcpy(&r, &pixelData[I], 4);
//		memcpy(&g, &pixelData[I+4], 4);
//		memcpy(&b, &pixelData[I+8], 4);

//		if(r==0.0 && g == 0.0 && b == 0.0) { r = 1.0; g = 1.0; b = 1.0; }

		glColor3f(pixelData[I]/255.F , pixelData[I+1]/255.F  , pixelData[I+2]/255.F );
		glVertex2f(x_coord, y_coord);

		// printf("%i %i %0.1f %0.1f %0.1f\n", x_coord, y_coord, r, g, b);

		x_coord++;
	}

	glEnd();

}

void GLUTCALLBACK Display (void) {
	glClear(GL_COLOR_BUFFER_BIT);
	ReceiveAndFillDispBuffers();
	glutSwapBuffers();
}




void OpenWindow (int pixels_x, int pixels_y) {

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(pixels_x, pixels_y);

	glutCreateWindow (" ");

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	glShadeModel(GL_FLAT);
	glDisable(GL_DITHER);

	glClear(GL_COLOR_BUFFER_BIT);

}



int main(int argc, char *argv[]) {

//	pixelData = (unsigned char *) malloc( sizeof(float) * 3 * pixels_x * pixels_y );
//	compressedData = (unsigned char *) malloc( sizeof(float) * 3 * pixels_x * pixels_y );

	pixelData = (unsigned char *) malloc( sizeof(unsigned char) * 3 * pixels_x * pixels_y );
	compressedData = (unsigned char *) malloc( sizeof(unsigned char) * 3 * pixels_x * pixels_y );

	struct hostent *he;
	struct sockaddr_in their_addr;

	if (argc != 2) {
		fprintf(stderr,"usage: client hostname\n");
		exit(1);
	}

	if ((he=gethostbyname(argv[1])) == NULL) {
		herror("gethostbyname");
		exit(1);
	}

	if ((sockfd = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
		perror("socket");
		exit(1);
	}

	their_addr.sin_family = AF_INET;
	their_addr.sin_port = htons(PORT);
	their_addr.sin_addr = *((struct in_addr *)he->h_addr);
	memset(their_addr.sin_zero, '\0', sizeof their_addr.sin_zero);

	if (connect(sockfd, (struct sockaddr *)&their_addr, sizeof their_addr) == -1) {
		perror("connect");
		exit(1);
	}

	glutInit(&argc, argv);
	OpenWindow(pixels_x, pixels_y);

	glLoadIdentity();
	gluOrtho2D(0.0, 512.0, 0.0, 512.0);
	glClearColor(1.0, 1.0, 1.0, 0.F);

	glutIdleFunc(Display);
	glutDisplayFunc(Display);
	glutMainLoop();

	close(sockfd);

	return 0;

} 

