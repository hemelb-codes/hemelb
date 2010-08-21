#ifndef HEME_VIS_VISTHREAD_H
#define HEME_VIS_VISTHREAD_H

#ifndef NO_STEER

#include <pthread.h>
#include <sys/types.h>

extern int doRendering;
extern int ShouldIRenderNow;
extern pthread_mutex_t var_lock;

extern u_int pixel_data_bytes;
extern u_int frame_details_bytes;
extern char* xdrSendBuffer_pixel_data;
extern char* xdrSendBuffer_frame_details;

void setRenderState(int val);
void* hemeLB_network(void *ptr);

#endif // NO_STEER

#endif // HEME_VIS_VISTHREAD_H
