#ifndef NO_STEER
#include <pthread.h>
#include <semaphore.h>

#include "rt.h"
#include "lb.h"

extern pthread_mutex_t network_buffer_copy_lock;
extern pthread_mutex_t LOCK;
extern pthread_cond_t network_send_frame;

extern sem_t nrl;
extern sem_t connected_sem;
extern sem_t steering_var_lock;

extern bool is_frame_ready;
extern bool connected;
extern bool sending_frame;

extern int send_array_length;


extern pthread_mutex_t steer_param_lock;


void  *hemeLB_steer (void*);
void UpdateSteerableParameters (int*, Vis*, LBM*);
#endif
