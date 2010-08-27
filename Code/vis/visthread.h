#ifndef HEMELB_VIS_VISTHREAD_H
#define HEMELB_VIS_VISTHREAD_H

#ifndef NO_STEER
#include <pthread.h>
#endif

#include <sys/types.h>

namespace hemelb
{
  namespace vis
  {
    extern int doRendering;
    extern int ShouldIRenderNow;
#ifndef NO_STEER
    extern pthread_mutex_t var_lock;
#endif

    extern u_int pixel_data_bytes;
    extern u_int frame_details_bytes;
    extern char* xdrSendBuffer_pixel_data;
    extern char* xdrSendBuffer_frame_details;

    void setRenderState(int val);
  }
}

#endif // HEMELB_VIS_VISTHREAD_H
