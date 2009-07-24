#ifndef NO_STEER


extern pthread_mutex_t steer_param_lock;


void  *hemeLB_steer (void*);
void UpdateSteerableParameters (int*, Vis*, LBM*);
#endif
