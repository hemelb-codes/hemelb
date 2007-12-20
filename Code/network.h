
//extern char *xdrSendBuffer_pixel_data;
//extern char *xdrSendBuffer_frame_details;

// data per pixel are colour id and pixel id (2 * sizeof(int) bytes)
extern int data_per_pixel;
extern int bytes_per_pixel_data;

// one int for colour_id and one for pixel id
extern u_int pixel_data_bytes;

// it is assumed that the frame size is the only detail
extern u_int frame_details_bytes;

extern int bits_per_char;
extern int bits_per_two_chars;

void setupNetworkBuffersAndThread();
void CleanUpNetworkBuffersAndThread();

void *hemeLB_network (void *ptr);

#ifdef RG

extern pthread_mutex_t network_buffer_copy_lock;
extern pthread_cond_t network_send_frame;

extern int send_array_length;

#endif // RG
