#ifndef NO_STEER

#include <pthread.h>

#include "config.h"
#include "visthread.h"

int doRendering = 0;
int ShouldIRenderNow = 0;

pthread_mutex_t var_lock = PTHREAD_MUTEX_INITIALIZER;

// data per pixel inlude the data for the pixel location and 4 colours
//in RGB format, thus #bytes per pixel are (sizeof(int)+4
//rgb)=sizeof(int)+4*3*sizeof(unsigned char))
int bytes_per_pixel_data = sizeof(int) + 4 * sizeof(unsigned char);

// one int for colour_id and one for pixel id
u_int pixel_data_bytes = COLOURED_PIXELS_MAX * bytes_per_pixel_data;

// it is assumed that the frame size is the only detail
u_int frame_details_bytes = 1 * sizeof(int);
char* xdrSendBuffer_pixel_data;
char* xdrSendBuffer_frame_details;

void setRenderState(int val)
{
  pthread_mutex_lock(&var_lock);
  //  doRendering = val;
  ShouldIRenderNow = val;
  pthread_mutex_unlock(&var_lock);
}

#endif // NO_STEER
