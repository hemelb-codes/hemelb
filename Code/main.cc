#ifndef NO_STEER

#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <sys/stat.h>

#include <sched.h>
#include <sys/param.h>

#include <semaphore.h>

#include "network.h"
#include "visthread.h"
#endif // NO_STEER

#include "config.h"
#include "steering.h"
#include "usage.h"
#include "benchmark.h"
#include "colourpalette.h"
#include "fileutils.h"
#include "utilityFunctions.h"

#define BCAST_FREQ   1

int cycle_id;
int time_step;
double intra_cycle_time;

#ifndef NO_STEER
bool updated_mouse_coords;
#endif

FILE *timings_ptr;

int main (int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output
  
#ifndef NO_STEER
  sem_init(&nrl, 0, 1);
  sem_init(&connected_sem, 0, 1);
  sem_init(&steering_var_lock, 0, 1);
  
  is_frame_ready = 0;
  connected = 0;
  sending_frame = 0;
  updated_mouse_coords = 0;
#endif
  
  double simulation_time;
  double minutes;
  double FS_time;
  double FS_plus_RT_time;
  double FS_plus_RT_plus_SL_time;
  
  int total_time_steps, stability = STABLE;
  int depths;
  int steering_session_id;
#ifdef NO_STEER
  int doRendering;
#endif
  
  int FS_time_steps;
  int FS_plus_RT_time_steps;
  int FS_plus_RT_plus_SL_time_steps;
  int snapshots_per_cycle, snapshots_period;
  int images_per_cycle, images_period;
  int is_unstable = 0;
  
#ifndef NO_STEER
  pthread_t network_thread;
  pthread_attr_t pthread_attrib;
#endif
  
  LBM lbm;
  
  Net net;
  
  SL sl;
  
#ifndef NOMPI
  int thread_level_provided;
  
  net.err = MPI_Init_thread (&argc, &argv, MPI_THREAD_FUNNELED, &thread_level_provided);
  net.err = MPI_Comm_size (MPI_COMM_WORLD, &net.procs);
  net.err = MPI_Comm_rank (MPI_COMM_WORLD, &net.id);
  
  if (net.id == 0)
    {
      printf("thread_level_provided %i\n", thread_level_provided);
    }
#else
  net.procs = 1;
  net.id = 0;
#endif
  
  check_conv = 0;
  
  if (argc == 5) // Check command line arguments
    {
      is_bench = 1;
      lbm.period     = atoi( argv[2] );
      lbm.voxel_size = atof( argv[3] );
      minutes        = atof( argv[4] );
    }
  else if (argc == 8)
    {
      is_bench = 0;
      lbm.cycles_max      = atoi( argv[2] );
      lbm.period          = atoi( argv[3] );
      lbm.voxel_size      = atof( argv[4] );
      snapshots_per_cycle = atoi( argv[5] );
      images_per_cycle    = atoi( argv[6] );
      steering_session_id = atoi( argv[7] );
      
      if (lbm.cycles_max > 1000)
	{
	  check_conv = 1;
	}
    }
  else
    {
      if (net.id == 0) Usage::printUsage(argv[0]);

#ifndef NOMPI
      net.err = MPI_Abort (MPI_COMM_WORLD, 1);
      net.err = MPI_Finalize ();
#else
      exit(1);
#endif
    }

  double total_time = UtilityFunctions::myClock();
  
  char* input_file_path( argv[1] );
  
  char input_config_name[256];
  char input_parameters_name[256];
  char vis_parameters_name[256];
  char timings_name[256];
  char procs_string[256];
  char snapshot_directory[256];
  char image_directory[256];
  char complete_image_name[256];
  
  strcpy ( input_config_name , input_file_path );
  strcat ( input_config_name , "/config.dat" );
  FileUtils::check_file(input_config_name);

  strcpy ( input_parameters_name , input_file_path );
  strcat ( input_parameters_name , "/pars.asc" );
  FileUtils::check_file(input_parameters_name);
  
  strcpy ( vis_parameters_name , input_file_path );
  strcat ( vis_parameters_name , "/rt_pars.asc" );
  FileUtils::check_file(vis_parameters_name);
  
  // Create directory path for the output images
  strcpy (image_directory, input_file_path);
  strcat (image_directory, "/Images/");
  
  //Create directory path for the output snapshots
  strcpy (snapshot_directory, input_file_path);
  strcat (snapshot_directory, "/Snapshots/");

  // Actually create the directories.
  if (net.id == 0)
  {
    FileUtils::MakeDirAllRXW(image_directory);
    FileUtils::MakeDirAllRXW(snapshot_directory);
  }
  
  sprintf ( procs_string, "%i", net.procs);
  strcpy ( timings_name , input_file_path );
  strcat ( timings_name , "/timings" );
  strcat ( timings_name , procs_string );
  strcat ( timings_name , ".asc" );

  if (net.id == 0)
    {
      timings_ptr = fopen (timings_name, "w");
    }
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "***********************************************************\n");
      fprintf (timings_ptr, "Opening parameters file:\n %s\n", input_parameters_name);
      fprintf (timings_ptr, "Opening config file:\n %s\n", input_config_name);
      fprintf (timings_ptr, "Opening vis parameters file:\n %s\n\n", vis_parameters_name);
    }
  
#ifndef NO_STEER
  if (!is_bench && net.id == 0)
    {
      xdrSendBuffer_pixel_data = (char *)malloc(pixel_data_bytes);
      xdrSendBuffer_frame_details = (char *)malloc(frame_details_bytes);

      pthread_mutex_init (&LOCK, NULL);
      pthread_cond_init (&network_send_frame, NULL);

  //    pthread_mutex_lock (&LOCK);
      
      pthread_attr_init (&pthread_attrib);
      pthread_attr_setdetachstate (&pthread_attrib, PTHREAD_CREATE_JOINABLE);
      
      pthread_create (&network_thread, &pthread_attrib, hemeLB_network, (void*)&steering_session_id);
    }
#endif // NO_STEER
  
  lbm.lbmInit (input_config_name, &net);
  
  lbm.lbmReadConfig (&net);
  
  lbm.lbmReadParameters (input_parameters_name, &net);
  
  if (netFindTopology (&net, &depths) == 0)
    {
      fprintf (timings_ptr, "MPI_Attr_get failed, aborting\n");
#ifndef NOMPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
  
  netInit (&lbm, &net);
  
  lbm.lbmSetInitialConditions (&net);
  
  visInit (&net, &vis, &sl);
  
  visReadParameters (vis_parameters_name, &lbm, &net, &vis);
#ifndef NO_STEER
  UpdateSteerableParameters (&doRendering, &vis, &lbm);
#endif

  FileUtils::DeleteDirContents (snapshot_directory);
  FileUtils::DeleteDirContents (image_directory);
  
  total_time_steps = 0;
  
  if (!is_bench)
    {
      int is_finished = 0;
      
      simulation_time = UtilityFunctions::myClock ();
      
      if (snapshots_per_cycle == 0)
	snapshots_period = 1e9;
      else
	snapshots_period = UtilityFunctions::max(1, lbm.period / snapshots_per_cycle);
      
      if (images_per_cycle == 0)
	images_period = 1e9;
      else
	images_period = UtilityFunctions::max(1, lbm.period / images_per_cycle);
      
      for (cycle_id = 1; cycle_id <= lbm.cycles_max && !is_finished; cycle_id++)
	{
	  int restart = 0;
	  
	  lbmInitMinMaxValues ();
	  
	  for (time_step = 1; time_step <= lbm.period; time_step++)
	    {
	      ++total_time_steps;
              intra_cycle_time = (PULSATILE_PERIOD * time_step) / lbm.period;
	      
	      int write_snapshot_image = (time_step % images_period == 0) ? 1 : 0;
#ifndef NO_STEER
	      int render_for_network_stream = 0;
	      
			/* In the following two if blocks we do the core magic to ensure we only render
			 when (1) we are not sending a frame or (2) we need to output to disk */
			
			if(net.id == 0) {		  
				sem_wait (&connected_sem);
				bool local_connected = connected;
				sem_post (&connected_sem);
				if(local_connected) {
					render_for_network_stream = (sending_frame == 0) ? 1 : 0;
				} else {
					render_for_network_stream = 0;
				}
			}

			if(total_time_steps%BCAST_FREQ == 0) {
				if(net.id == 0)
					doRendering = (render_for_network_stream || write_snapshot_image) ? 1 : 0;
				if(net.id == 0)
					sem_wait (&steering_var_lock);
				UpdateSteerableParameters (&doRendering, &vis, &lbm);
				if(net.id == 0)
					sem_post (&steering_var_lock);
			}

			/* for debugging purposes we want to ensure we capture the variables in a single
			 instant of time since variables might be altered by the thread half way through?
			 This is to be done. */			
			
			if(net.id == 0 && time_step%100==0)
				printf("time step %i sending_frame %i render_network_stream %i write_snapshot_image %i rendering %i\n",
					   time_step, sending_frame, render_for_network_stream, write_snapshot_image, doRendering);

#endif // NO_STEER

	      lbmUpdateBoundaryDensities (cycle_id, time_step, &lbm);
	      
	      if (!check_conv)
		{
		  stability = lbmCycle (doRendering, &lbm, &net);
		  
		  if ((restart = lbmIsUnstable (&net)) != 0)
		    {
		      break;
		    }
		  lbmUpdateInletVelocities (time_step, &lbm, &net);
		}
	      else
		{
		  stability = lbmCycle (cycle_id, time_step, doRendering, &lbm, &net);
		  
		  if (stability == UNSTABLE)
		    {
		      restart = 1;
		      break;
		    }
		  lbmUpdateInletVelocities (time_step, &lbm, &net);
		}
#ifndef NO_STREAKLINES
	      slStreakLines (time_step, lbm.period, &net, &sl);
#endif
#ifndef NO_STEER
	      if (total_time_steps%BCAST_FREQ == 0 && doRendering && !write_snapshot_image)
		{
		  visRender (RECV_BUFFER_A, ColourPalette::PickColour, &net, &sl);
		  
		  if (vis_mouse_x >= 0 && vis_mouse_y >= 0 && updated_mouse_coords)
		    {
		      for (int i = 0; i < col_pixels_recv[RECV_BUFFER_A]; i++)
			{
			  if ((col_pixel_recv[RECV_BUFFER_A][i].i & RT) &&
			      (col_pixel_recv[RECV_BUFFER_A][i].i & PIXEL_ID_MASK) == PixelId (vis_mouse_x,vis_mouse_y))
			    {
			      visCalculateMouseFlowField (&col_pixel_recv[RECV_BUFFER_A][i], &lbm);
			      break;
			    }
			}
		      updated_mouse_coords = 0;
		    }
		 if (net.id == 0)
		   {
		     is_frame_ready = 1;
		     sem_post(&nrl); // let go of the lock
		   }
		}
#endif // NO_STEER
	      if (write_snapshot_image)
		{
		  visRender (RECV_BUFFER_B, ColourPalette::PickColour, &net, &sl);
		  
		  if (net.id == 0)
		    {
		      char image_filename[255];
		      
		      snprintf(image_filename, 255, "%08i.dat", time_step);
		      strcpy ( complete_image_name, image_directory );
		      strcat ( complete_image_name, image_filename );
		      
		      visWriteImage (RECV_BUFFER_B, complete_image_name, ColourPalette::PickColour);
		    }
		}
	      if (time_step%snapshots_period == 0)
		{
		  char snapshot_filename[255];
		  char complete_snapshot_name[255];
		  
		  snprintf(snapshot_filename, 255, "snapshot_%06i.asc", time_step);
		  strcpy ( complete_snapshot_name, snapshot_directory );
		  strcat ( complete_snapshot_name, snapshot_filename );
		  
		  lbmWriteConfigASCII (stability, complete_snapshot_name, &lbm, &net);
		}
#ifndef NO_STEER
	      if (net.id == 0)
		{
                  if (render_for_network_stream == 1)
		    {
		      // printf("sending signal to thread that frame is ready to go...\n"); fflush(0x0);
		      sched_yield();
		      sem_post( &nrl );
		      //pthread_mutex_unlock (&LOCK);
		      //pthread_cond_signal (&network_send_frame);
		    }
		}
#endif
	      if (stability == STABLE_AND_CONVERGED)
		{
		  is_finished = 1;
		  break;
		}
	      if (lbm_terminate_simulation)
		{
		  is_finished = 1;
		  break;
		}
	      if (lbm.period > 400000)
		{
		  is_unstable = 1;
		  break;
		}
	    }
	  
	  if (restart)
	    {
	      FileUtils::DeleteDirContents (snapshot_directory);
              FileUtils::DeleteDirContents (image_directory);
	      
	      lbm.lbmRestart (&net);
#ifndef NO_STREAKLINES
	      slRestart (&sl);
#endif
	      if (net.id == 0)
		{
		  printf ("restarting: period: %i\n", lbm.period);
		  fflush (0x0);
		}
	      snapshots_period = (snapshots_per_cycle == 0) ? 1e9 : UtilityFunctions::max(1, lbm.period/snapshots_per_cycle);
	      
	      images_period = (images_per_cycle == 0) ? 1e9 : UtilityFunctions::max(1, lbm.period/images_per_cycle);
	      
	      cycle_id = 0;
	      continue;
	    }
	  lbmCalculateFlowFieldValues (&lbm);
	  
	  if (net.id == 0)
	    {
	      if (!check_conv)
		{
		  fprintf (timings_ptr, "cycle id: %i\n", cycle_id);
		  printf ("cycle id: %i\n", cycle_id);
		}
	      else
		{
		  fprintf (timings_ptr, "cycle id: %i, conv_error: %le\n", cycle_id, conv_error);
		  printf ("cycle id: %i, conv_error: %le\n", cycle_id, conv_error);
		}
	      fflush(NULL);
	    }
	}
      simulation_time = UtilityFunctions::myClock () - simulation_time;
      
      time_step = UtilityFunctions::min(time_step, lbm.period);
      cycle_id = UtilityFunctions::min(cycle_id, lbm.cycles_max);
      time_step = time_step * cycle_id;
    }
  else // is_bench
    {
      double elapsed_time;
      
      int bench_period = (int)fmax(1.0, (1e+6 * net.procs) / lbm.total_fluid_sites);
      
      // benchmarking HemeLB's fluid solver only
      
      FS_time = UtilityFunctions::myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  ++total_time_steps;
	  lbmUpdateBoundaryDensities (total_time_steps/lbm.period, total_time_steps%lbm.period, &lbm);
	  stability = lbmCycle (0, &lbm, &net);
	  
	  elapsed_time = UtilityFunctions::myClock () - FS_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " FS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      BenchmarkTimer::IsBenchSectionFinished (1.0, elapsed_time))
	    {
	      break;
	    }
	}

      FS_time_steps = (int)(time_step * minutes / (3 * 1.0) - time_step);
      FS_time = UtilityFunctions::myClock ();
      
      for (time_step = 1; time_step <= FS_time_steps; time_step++)
	{
	  ++total_time_steps;
	  lbmUpdateBoundaryDensities (total_time_steps/lbm.period, total_time_steps%lbm.period, &lbm);
	  stability = lbmCycle (1, &lbm, &net);
	}
      FS_time = UtilityFunctions::myClock () - FS_time;
      
      
      // benchmarking HemeLB's fluid solver and ray tracer
      
      vis_mode = 0;
      vis_image_freq = 1;
      vis_streaklines = 0;
      FS_plus_RT_time = UtilityFunctions::myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  ++total_time_steps;
	  lbmUpdateBoundaryDensities (total_time_steps/lbm.period, total_time_steps%lbm.period, &lbm);
	  stability = lbmCycle (1, &lbm, &net);
	  visRender (RECV_BUFFER_A, ColourPalette::PickColour, &net, &sl);
	  
	  // partial timings
	  elapsed_time = UtilityFunctions::myClock () - FS_plus_RT_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " FS + RT, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      BenchmarkTimer::IsBenchSectionFinished (1.0, elapsed_time))
	    {
	      break;
	    }
	}
      FS_plus_RT_time_steps = (int)(time_step * minutes / (3 * 1.0) - time_step);
      FS_plus_RT_time = UtilityFunctions::myClock ();
      
      for (time_step = 1; time_step <= FS_plus_RT_time_steps; time_step++)
	{
	  ++total_time_steps;
	  lbmUpdateBoundaryDensities (total_time_steps/lbm.period, total_time_steps%lbm.period, &lbm);
	  stability = lbmCycle (1, &lbm, &net);
	  visRender (RECV_BUFFER_A, ColourPalette::PickColour, &net, &sl);
	}
      FS_plus_RT_time = UtilityFunctions::myClock () - FS_plus_RT_time;
      
#ifndef NO_STREAKLINES
      // benchmarking HemeLB's fluid solver, ray tracer and streaklines
      
      vis_mode = 2;
      vis_streaklines = 1;
      FS_plus_RT_plus_SL_time = UtilityFunctions::myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  ++total_time_steps;
	  lbmUpdateBoundaryDensities (total_time_steps/lbm.period, total_time_steps%lbm.period, &lbm);
	  stability = lbmCycle (1, &lbm, &net);
          slStreakLines (time_step, lbm.period, &net, &sl);
	  visRender (RECV_BUFFER_A, ColourPalette::PickColour, &net, &sl);
	  
	  // partial timings
	  elapsed_time = UtilityFunctions::myClock () - FS_plus_RT_plus_SL_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " FS + RT + SL, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      BenchmarkTimer::IsBenchSectionFinished (1.0, elapsed_time))
	    {
	      break;
	    }
	}
      FS_plus_RT_plus_SL_time_steps = (int)(time_step * minutes / (3 * 1.0) - time_step);
      FS_plus_RT_plus_SL_time = UtilityFunctions::myClock ();
      
      for (time_step = 1; time_step <= FS_plus_RT_plus_SL_time_steps; time_step++)
	{
	  ++total_time_steps;
	  lbmUpdateBoundaryDensities (total_time_steps/lbm.period, total_time_steps%lbm.period, &lbm);
	  stability = lbmCycle (1, &lbm, &net);
	  slStreakLines (time_step, lbm.period, &net, &sl);
	  visRender (RECV_BUFFER_A, ColourPalette::PickColour, &net, &sl);
	}
      FS_plus_RT_plus_SL_time = UtilityFunctions::myClock () - FS_plus_RT_plus_SL_time;
#endif // NO_STREAKLINES
    } // is_bench
  
  if (!is_bench)
    {  
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "\n");
	  fprintf (timings_ptr, "threads: %i, machines checked: %i\n\n", net.procs, net_machines);
	  fprintf (timings_ptr, "topology depths checked: %i\n\n", depths);
	  fprintf (timings_ptr, "fluid sites: %i\n\n", lbm.total_fluid_sites);
	  fprintf (timings_ptr, "cycles and total time steps: %i, %i \n\n", cycle_id, total_time_steps);
	  fprintf (timings_ptr, "time steps per second: %.3f\n\n", total_time_steps / simulation_time);
	}
    }
  else  // is_bench
    {
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "\n---------- BENCHMARK RESULTS ----------\n");
	  
	  fprintf (timings_ptr, "threads: %i, machines checked: %i\n\n", net.procs, net_machines);
	  fprintf (timings_ptr, "topology depths checked: %i\n\n", depths);
	  fprintf (timings_ptr, "fluid sites: %i\n\n", lbm.total_fluid_sites);
	  fprintf (timings_ptr, " FS, time steps per second: %.3f, MSUPS: %.3f, time: %.3f\n\n",
		   FS_time_steps / FS_time,
		   1.e-6 * lbm.total_fluid_sites / (FS_time / FS_time_steps),
		   FS_time);
	  
	  fprintf (timings_ptr, " FS + RT, time steps per second: %.3f, time: %.3f\n\n",
		   FS_plus_RT_time_steps / FS_plus_RT_time, FS_plus_RT_time);
#ifndef NO_STREAKLINES
	  fprintf (timings_ptr, " FS + RT + SL, time steps per second: %.3f, time: %.3f\n\n",
		   FS_plus_RT_plus_SL_time_steps / FS_plus_RT_plus_SL_time, FS_plus_RT_plus_SL_time);
#endif
	}
    }
  
  if (is_unstable)
    {
      if (net.id == 0)
	{
	  fprintf (timings_ptr, "Attention: simulation unstable with %i timesteps/cycle\n",
		   lbm.period);
	  fprintf (timings_ptr, "Simulation is terminated\n");
	  fclose (timings_ptr);
	}
    }
  else
    {
      if (net.id == 0)
	{
	  if (!is_bench)
	    {
	      vis_pressure_min = lbm.lbmConvertPressureToPhysicalUnits (lbm_density_min * Cs2);
	      vis_pressure_max = lbm.lbmConvertPressureToPhysicalUnits (lbm_density_max * Cs2);
	      
	      vis_velocity_min = lbm.lbmConvertVelocityToPhysicalUnits (lbm_velocity_min);
	      vis_velocity_max = lbm.lbmConvertVelocityToPhysicalUnits (lbm_velocity_max);
	      
	      vis_stress_min = lbm.lbmConvertStressToPhysicalUnits (lbm_stress_min);
	      vis_stress_max = lbm.lbmConvertStressToPhysicalUnits (lbm_stress_max);
	      
	      fprintf (timings_ptr, "time steps per cycle: %i\n", lbm.period);
	      fprintf (timings_ptr, "pressure min, max (mmHg): %le, %le\n", vis_pressure_min, vis_pressure_max);
	      fprintf (timings_ptr, "velocity min, max (m/s) : %le, %le\n", vis_velocity_min, vis_velocity_max);
	      fprintf (timings_ptr, "stress   min, max (Pa)  : %le, %le\n", vis_stress_min, vis_stress_max);
	      fprintf (timings_ptr, "\n");
	      
	      for (int n = 0; n < lbm.inlets; n++)
		{
		  fprintf (timings_ptr, "inlet id: %i, average / peak velocity (m/s): %le / %le\n",
			   n, lbm_average_inlet_velocity[ n ], lbm_peak_inlet_velocity[ n ]);
		}
	      fprintf (timings_ptr, "\n");
	    }
	  fprintf (timings_ptr, "\n");
	  fprintf (timings_ptr, "domain decomposition time (s):             %.3f\n", net.dd_time);
	  fprintf (timings_ptr, "pre-processing buffer management time (s): %.3f\n", net.bm_time);
	  fprintf (timings_ptr, "input configuration reading time (s):      %.3f\n", net.fr_time);
	  
	  total_time = UtilityFunctions::myClock () - total_time;
	  fprintf (timings_ptr, "total time (s):                            %.3f\n\n", total_time);
	  
	  fprintf (timings_ptr, "Sub-domains info:\n\n");
	  
	  for (int n = 0; n < net.procs; n++)
	    {
	      fprintf (timings_ptr, "rank: %i, fluid sites: %i\n", n, net.fluid_sites[ n ]);
	    }
	  
	  fclose (timings_ptr);
	}
    }
  visEnd (&sl);
  netEnd (&net);
  lbm.lbmEnd ();
  
#ifndef NO_STEER
  if (!is_bench && net.id == 0)
    {
      // there are some problems if the following function is called
      
      //pthread_join (network_thread, NULL);
      free(xdrSendBuffer_frame_details);
      free(xdrSendBuffer_pixel_data);
    }
#endif
  
#ifndef NOMPI
  net.err = MPI_Finalize ();
#endif
  
  return(0);
}

