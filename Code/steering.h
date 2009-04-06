#ifndef NO_STEER

extern pthread_mutex_t steer_param_lock;
void  *hemeLB_steer (void*);

#endif


void UpdateSteerableParameters (int*, Vis*, LBM*);

