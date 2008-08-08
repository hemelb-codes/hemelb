#include <semaphore.h>
#include <pthread.h>
#include <unistd.h>

#include "config.h"
#include "steering.h"
#include "network.h"
#include "visthread.h"
#include "colourpalette.h"

// RealityGrid includes
#include "ReG_Steer_types.h"
#include "ReG_Steer_Appside.h"

pthread_mutex_t reg_lock = PTHREAD_MUTEX_INITIALIZER;

// steering parameters
char* steer_par_label[STEERABLE_PARAMETERS-1] = {"ctr_x", "ctr_y", "ctr_z",
						 "longitude", "latitude",
						 "zoom", "vis_brightness",
						 "velocity_max", "stress_max",
						 "mouse_x", "mouse_y"};
// extra params here as we need to initialize the "terminate"
// and "render" params even if we don't them.
float steer_par[STEERABLE_PARAMETERS+1] = {0.0, 0.0, 0.0, 45.0, 45.0, 1.0,
					 0.03, 0.1, 0.1, -1.0, -1.0, 0.0, 1.0};
char* steer_par_min[STEERABLE_PARAMETERS-1] = {"0.0", "0.0", "0.0", "-90.0",
					       "-90.0", "0.5", "", "", "",
					       "", ""};
char* steer_par_max[STEERABLE_PARAMETERS-1] = {"1.0", "1.0", "1.0", "90.0",
					       "90.0", "5.0", "", "", "",
					       "", ""};
int   reg_finished;
sem_t reg_wait;

void* hemeLB_network(void *ptr) {

  setRenderState(0);

  printf("kicking off RealityGrid network thread.....\n");
  fflush(0x0);

  int reg_status;
  int reg_num_cmds;
  int reg_cmds[REG_INITIAL_NUM_CMDS];
  int reg_io_type;
  int reg_io_channel;

  pthread_t steering_thread;
  pthread_attr_t steering_thread_attrib; 
  pthread_attr_init(&steering_thread_attrib);
  pthread_attr_setdetachstate(&steering_thread_attrib,
			      PTHREAD_CREATE_JOINABLE);

  // initialize done signals and semaphore
  reg_finished = 0;
  sem_init(&reg_wait, 0, 0);
  int frame_number = 0;
  int done = 0;
  //unsigned int pix_data[4];
  unsigned int* pix_index = NULL;
  unsigned char* pix_colours = NULL;
  int old_col_pixels = col_pixels;

  pthread_mutex_lock(&LOCK);

  // initialize steering library
  Steering_enable(REG_TRUE);
      
  reg_num_cmds = 2;
  reg_cmds[0] = REG_STR_STOP;
  reg_cmds[1] = REG_STR_PAUSE_INTERNAL;

  reg_status = Steering_initialize ("HemeLB", reg_num_cmds, reg_cmds);

  // TODO
  // what do we *really* do on failure here?
  if(reg_status == REG_FAILURE) {
    printf("Failed to initialize RealityGrid steering\n");
    return 0;
  }

  // register parameters
  for(int i = 0; i < (STEERABLE_PARAMETERS-1); i++) {
    reg_status = Register_param(steer_par_label[i], REG_TRUE,
				(void*)(&steer_par[i]), REG_FLOAT,
				steer_par_min[i], steer_par_max[i]);

    if(reg_status != REG_SUCCESS) break;
  }

  // TODO
  // what do we *really* do on failure here?
  if(reg_status != REG_SUCCESS) {
    printf("Failed to register RealityGrid steering parameters\n");
    return 0;
  }

  // register IO channel
  reg_status = Register_IOType("HemeLB RT Viz Data",
			       REG_IO_OUT, 1, &reg_io_type);

  // TODO
  // what do we *really* do on failure here?
  if(reg_status != REG_SUCCESS) {
    printf("Failed to register RealityGrid IO channel\n");
    return 0;
  }

  // create steering thread
  pthread_create(&steering_thread, &steering_thread_attrib,
		 hemeLB_steer, (void*) &reg_io_type);


  pthread_mutex_unlock(&LOCK);

  setRenderState(1);


  while(!done) {

    printf("THREAD: waiting for signal that frame is ready to send..\n");
    fflush(0x0);

    pthread_mutex_lock(&LOCK);
    pthread_cond_wait(&network_send_frame, &LOCK);
    setRenderState(0);

    printf("THREAD: received signal that frame is ready to send..\n");
    fflush(0x0);

    // try to send frame down ReG pipe
    if(Emit_start(reg_io_type, frame_number, &reg_io_channel) == REG_SUCCESS) {
      if(old_col_pixels != col_pixels) {
	pix_index = (unsigned int*) realloc(pix_index, sizeof(unsigned int) * col_pixels);
	pix_colours = (unsigned char*) realloc(pix_colours, sizeof(unsigned char) * col_pixels * 12);
	old_col_pixels = col_pixels;
      }

      for(int i = 0; i < col_pixels; i++) {
	rawWritePixel(&col_pixel_recv[i], &pix_index[i],
		      &pix_colours[(i * 12)], ColourPalette);
      }

      
      reg_status = Emit_data_slice(reg_io_channel, REG_INT,
    				   col_pixels, (void*) pix_index);
      reg_status = Emit_data_slice(reg_io_channel, REG_CHAR,
    				   (col_pixels * 12), (void*) pix_colours);
      Emit_stop(&reg_io_channel);
    }

    setRenderState(1);
    pthread_mutex_unlock(&LOCK);  
    frame_number++;

    // finished?
    pthread_mutex_lock(&reg_lock);
    done = reg_finished;
    pthread_mutex_unlock(&reg_lock);
  } // !done

  // tell steering thread that this one is done
  sem_post(&reg_wait);
}

void* hemeLB_steer(void* ptr) {
  int    reg_status;
  int    reg_num_params_changed;
  int    reg_num_recvd_cmds;
  int    reg_recvd_cmds[REG_MAX_NUM_STR_CMDS];
  char** reg_changed_param_labels;
  char** reg_recvd_cmd_params;  
  int    reg_timestep;
  int    reg_io_type = *((int*)ptr);

  reg_changed_param_labels = Alloc_string_array(REG_MAX_STRING_LENGTH,
						REG_MAX_NUM_STR_PARAMS);
  reg_recvd_cmd_params = Alloc_string_array(REG_MAX_STRING_LENGTH,
					    REG_MAX_NUM_STR_CMDS);

  printf("Kicking off RealityGrid steering thread\n");

  // steering control loop
  reg_timestep = 0;

  while(steer_par[STEERABLE_PARAMETERS-1] == 0.0) {
    usleep(1000);
    printf("REG: Steer loop %d\n", reg_timestep);
    fflush(0);
    reg_status = Steering_control(reg_timestep++,
				  &reg_num_params_changed,
				  reg_changed_param_labels,
				  &reg_num_recvd_cmds,
				  reg_recvd_cmds,
				  reg_recvd_cmd_params);

    // ignore errors here and carry on to the next loop
    if(reg_status != REG_SUCCESS) {
      printf("Steering_control failed\n");
      continue;
    }

    // process commands received
    for(int i = 0; i < reg_num_recvd_cmds; i++) {
      switch(reg_recvd_cmds[i]) {
      case REG_STR_STOP:
	printf("Received STOP command\n");

	steer_par[STEERABLE_PARAMETERS-1] = 1.0;
	pthread_mutex_lock(&reg_lock);
	reg_finished = 1;
	pthread_mutex_unlock(&reg_lock);

	break;

      default:
	if(reg_recvd_cmds[i] == reg_io_type) {
	  printf("output vis!\n");
	  fflush(0);
	}
	break;
      }
    } // end of command processing

  } // end of steering loop

  // wait on semaphore until vis output loop has finished
  printf("REG: Waiting to clean up...\n");
  fflush(0);
  sem_wait(&reg_wait);

  // clean up steering
  Steering_finalize();
  printf("REG: Cleaned up...\n");
  fflush(0);

  return 0;
}
