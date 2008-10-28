#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#include <sys/stat.h>

#include "config.h"
#include "network.h"
#include "steering.h"
#include "usage.h"
#include "benchmark.h"
#include "colourpalette.h"
#include "visthread.h"
#include "fileutils.h"

FILE* timings_ptr;

int main (int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output
  
  double simulation_time;
  double minutes;
  double fluid_solver_time;
  double fluid_solver_and_vis_time;
  double vis_without_compositing_time;
  double steering_and_vis_time;
  double io_time, other_time;
  double start_time, end_time;
  
  int cycle_id;
  int total_time_steps, time_step, stability = STABLE;
  int depths;
  
  int fluid_solver_time_steps;
  int fluid_solver_and_vis_time_steps;
  int vis_without_compositing_time_steps;
  int snapshots_per_cycle, snapshots_period;
  int images_per_cycle, images_period;
  
  pthread_t network_thread;
  pthread_attr_t pthread_attrib;
  
  LBM lbm;
  
  Net net;
  
  SL sl;
  
#ifndef NOMPI
  net.err = MPI_Init (&argc, &argv);
  net.err = MPI_Comm_size (MPI_COMM_WORLD, &net.procs);
  net.err = MPI_Comm_rank (MPI_COMM_WORLD, &net.id);
#else
  net.procs = 1;
  net.id = 0;
#endif
  
  check_conv = 0;
  
  if (argc == 3) // Check command line arguments
    {
      is_bench = 1;
      minutes = atof( argv[2] );
    }
  else if (argc == 7)
    {
      is_bench = 0;
      lbm.cycles_max      = atoi( argv[2] );
      lbm.period          = atoi( argv[3] );
      lbm.voxel_size      = atof( argv[4] );
      snapshots_per_cycle = atoi( argv[5] );
      images_per_cycle    = atoi( argv[6] );
      
      if (lbm.cycles_max > 1000)
	{
	  check_conv = 1;
	}
    }
  else
    {
      if (net.id == 0) usage(argv[0]);

#ifndef NOMPI
      net.err = MPI_Abort (MPI_COMM_WORLD, 1);
      net.err = MPI_Finalize ();
#else
      exit(1);
#endif
    }

  double total_time = myClock();
  
  char* input_file_path( argv[1] );
  
  char input_config_name[256];
  char input_parameters_name[256];
  char output_config_name[256];
  char vis_parameters_name[256];
  char timings_name[256];
  char procs_string[256];
  char snapshot_directory[256];
  char image_directory[256];
  char complete_image_name[256];
  char rm_files[256];
  
  strcpy ( input_config_name , input_file_path );
  strcat ( input_config_name , "/config.dat" );
  check_file(input_config_name);

  strcpy ( input_parameters_name , input_file_path );
  strcat ( input_parameters_name , "/pars.asc" );
  check_file(input_parameters_name);
  
  strcpy ( vis_parameters_name , input_file_path );
  strcat ( vis_parameters_name , "/rt_pars.asc" );
  check_file(vis_parameters_name);

  strcpy ( output_config_name , input_file_path );
  strcat ( output_config_name , "/out.dat" );
  
  // Create directory for the output images
  strcpy (image_directory, input_file_path);
  strcat (image_directory, "/Images/");
  if (net.id == 0) mkdir  (image_directory, 0777);
  
  //Create directory for the output snapshots
  strcpy (snapshot_directory, input_file_path);
  strcat (snapshot_directory, "/Snapshots/");
  if (net.id == 0) mkdir  (snapshot_directory, 0777);
  
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
  
  if(net.id == 0)
    {
      xdrSendBuffer_pixel_data = (char *)malloc(pixel_data_bytes);
      xdrSendBuffer_frame_details = (char *)malloc(frame_details_bytes);

      pthread_mutex_init (&LOCK, NULL);
      pthread_cond_init (&network_send_frame, NULL);

  //    pthread_mutex_lock (&LOCK);
      
      pthread_attr_init (&pthread_attrib);
      pthread_attr_setdetachstate (&pthread_attrib, PTHREAD_CREATE_JOINABLE);
      
      pthread_create (&network_thread, &pthread_attrib, hemeLB_network, NULL);
    }

  lbmReadParameters (input_parameters_name, &lbm, &net);

  lbmInit (input_config_name, &lbm, &net);
  
  if (netFindTopology (&net, &depths) == 0)
    {
      fprintf (timings_ptr, "MPI_Attr_get failed, aborting\n");
#ifndef NOMPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
  
  netInit (&lbm, &net);
  
  lbmSetInitialConditions (&net);
  
  visInit (&net, &vis, &sl);
  
  visReadParameters (vis_parameters_name, &lbm, &net, &vis);
  
  // OUCH!!!! This needs to be changed - SM
  if (net.id == 0)
    {
      strcpy (rm_files, "rm ");
      strcat (rm_files, snapshot_directory);
      strcat (rm_files, "*");
      system (rm_files);
      
      strcpy (rm_files, "rm ");
      strcat (rm_files, image_directory);
      strcat (rm_files, "*");
      system (rm_files);
    }
  
  if (!is_bench)
    {
      int is_finished = 0;
      
      total_time_steps = 0;
      
      simulation_time = myClock ();
      fluid_solver_time = 0.;
      steering_and_vis_time = myClock ();
      io_time = 0.0;
      other_time = 0.0;
   
      if (snapshots_per_cycle == 0)
	snapshots_period = 1e9;
      else
	snapshots_period = max(1, lbm.period / snapshots_per_cycle);
      
      if (images_per_cycle == 0)
	images_period = 1e9;
      else
	images_period = max(1, lbm.period / images_per_cycle);
      
      for (cycle_id = 1; cycle_id <= lbm.cycles_max && !is_finished; cycle_id++)
	{
	  int restart = 0;
	  
	  for (time_step = 1; time_step <= lbm.period; time_step++, total_time_steps++)
	    {

              int render_for_network_stream = 0;
	      int write_snapshot_image = 0;

	      if (time_step % images_period == 0) {
                write_snapshot_image = 1;
                doRendering = 1;
		printf("Setting write_snapshot_image -> 1, time step %i\n", time_step);
              }

	      if (net.id == 0)
		{
		  int lock_return = pthread_mutex_trylock ( &LOCK );
		  
		  // printf("attempting to aquire mutex lock -> %i -> ", lock_return); 
		  
		  if (lock_return == EBUSY) {
		    // printf("lock busy\n");
		  } else {
		    // printf("aquired lock\n");
		    render_for_network_stream = 1;
		printf("Setting render_for_network_stream -> 1, time step %i\n", time_step);
		    doRendering = 1;
		  }
		  // printf("ShouldIRenderNow %i\n", ShouldIRenderNow); fflush(0x0);
		}

	      UpdateSteerableParameters (&doRendering, &vis, &lbm);
	      
	      // printf (" doRendering: %i\n", doRendering); Between
	      // the visRenderA/B calls, do not change any vis
	      // parameters.
	      
	      // if(setRendering==1) { doRendering = 1; setRendering = 0; }

	      if (doRendering)
		{
		  visRenderA (ColourPalette, &net, &sl);
		}

	      lbmVaryBoundaryDensities (cycle_id, time_step, &lbm);
	      
	      if (!check_conv)
		{
		  start_time = myClock ();
		  
		  stability = lbmCycle (cycle_id, time_step, doRendering, &lbm, &net);
		  
		  if ((restart = lbmIsUnstable (&net)) != 0)
		    {
		      end_time = myClock ();
		      fluid_solver_time += end_time - start_time;
		      break;
		    }
		  lbmUpdateInletFluxes (time_step, &lbm, &net);
		  
		  end_time = myClock ();
		  fluid_solver_time += end_time - start_time;
		}
	      else
		{
		  start_time = myClock ();
		  
		  stability = lbmCycleConv (cycle_id, time_step, doRendering, &lbm, &net);
		  
		  end_time = myClock ();
		  fluid_solver_time += end_time - start_time;
		  
		  if (stability == UNSTABLE)
		    {
		      restart = 1;
		      break;
		    }
		  lbmUpdateInletFluxes (time_step, &lbm, &net);
		  
		  end_time = myClock ();
		  fluid_solver_time += end_time - start_time;
		}

	      if (write_snapshot_image)
		{
		  start_time = myClock ();
		  
		  char image_filename[255];
		  
		  snprintf(image_filename, 255, "%08i.dat", time_step);
		  strcpy ( complete_image_name, image_directory );
		  strcat ( complete_image_name, image_filename );
		  
		  end_time = myClock ();
		  io_time += end_time - start_time;
		}

	      slStreakLines (time_step, lbm.period, &net, &sl);
	      
	      if (doRendering)
		{
		  visRenderB (write_snapshot_image, complete_image_name, ColourPalette, &net);
		  
		  for (int i = 0; i < col_pixels; i++)
		    {
		      if ((col_pixel_recv[ i ].i & RT) &&
			  (col_pixel_recv[ i ].i & PIXEL_ID_MASK) == PixelId (vis_mouse_x, vis_mouse_y))
			{
			  visCalculateMouseFlowField (&col_pixel_recv[ i ], &lbm);
			}
		    }
		}

	      if (time_step%snapshots_period == 0)
		{
		  start_time = myClock ();
		  
		  char snapshot_filename[255];
		  char complete_snapshot_name[255];
		  
		  snprintf(snapshot_filename, 255, "snapshot_%06i.asc", time_step);
		  strcpy ( complete_snapshot_name, snapshot_directory );
		  strcat ( complete_snapshot_name, snapshot_filename );
		  
		  lbmWriteConfigASCII (stability, complete_snapshot_name, &lbm, &net);
		  
		  end_time = myClock ();
		  io_time += end_time - start_time;
		}

	      if (net.id == 0)
		{
                  if( render_for_network_stream == 1 ) {
		    // printf("sending signal to thread that frame is ready to go...\n"); fflush(0x0);
		    pthread_mutex_unlock (&LOCK);
		    pthread_cond_signal (&network_send_frame);
		    // printf("...signal sent\n"); fflush(0x0);
                  }
                  if( doRendering == 1 ) doRendering = 0;
		} 

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
	      if (net.id == 0)
		{
		  if (time_step%1 == 0)
		    printf ("time step: %i\n", time_step);
		}
	    }

	  if (restart)
	    {
	      start_time = myClock ();
	      
	      lbmRestart (&lbm, &net);
	      
	      slRestart (&sl);
	      
	      if (net.id == 0)
		{
		  strcpy (rm_files, "rm ");
		  strcat (rm_files, snapshot_directory);
		  strcat (rm_files, "*");
		  system (rm_files);
		  
		  strcpy (rm_files, "rm ");
		  strcat (rm_files, image_directory);
		  strcat (rm_files, "*");
		  system (rm_files);
		}
	      if (net.id == 0)
		{
		  printf ("restarting: period: %i\n", lbm.period);
		  fflush (0x0);
		}
	      if (snapshots_per_cycle == 0)
		snapshots_period = 1e9;
	      else
		snapshots_period = max(1, lbm.period / snapshots_per_cycle);
	      
	      if (images_per_cycle == 0)
		images_period = 1e9;
	      else
		images_period = max(1, lbm.period / images_per_cycle);
	      
	      cycle_id = 1;
	      
	      end_time = myClock ();
	      other_time += (end_time - start_time);
	      continue;
	    }

	  start_time = myClock ();
	  
	  lbmCalculateFlowFieldValues (cycle_id, lbm.period, &lbm);
	  
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
	  end_time = myClock ();
	  other_time += end_time - start_time;
	}
      
      if (net.procs > 1)
	{
	  if (net.id == 1)
	    {
	      MPI_Send (&fluid_solver_time, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
	    }
	  else if (net.id == 0)
	    {
	      MPI_Recv (&fluid_solver_time, 1, MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, net.status);
	    }
	}

      steering_and_vis_time = myClock () - steering_and_vis_time - fluid_solver_time - io_time-other_time;
      simulation_time = myClock () - simulation_time;
      time_step = (min(time_step, lbm.period)) * min(cycle_id, lbm.cycles_max);

    }
  else // is_bench
    {
      double elapsed_time;
  
      int bench_period = (int)fmax(1., (1e+6 * net.procs) / lbm.total_fluid_sites);
      
      // benchmarking HemeLB's fluid solver only
      
      fluid_solver_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  stability = lbmCycle (1, 1, 0, &lbm, &net);
	  
	  // partial timings
	  elapsed_time = myClock () - fluid_solver_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " FS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      IsBenchSectionFinished (0.5, elapsed_time))
	    {
	      break;
	    }
	}

      fluid_solver_time_steps = (int)(time_step * minutes / (3. * 0.5) - time_step);
      fluid_solver_time = myClock ();
      
      for (time_step = 1; time_step <= fluid_solver_time_steps; time_step++)
	{
	  stability = lbmCycle (1, 1, 1, &lbm, &net);
	}

      fluid_solver_time = myClock () - fluid_solver_time;
      
      
      // benchmarking HemeLB's fluid solver and ray tracer
      
      vis_image_freq = 1;
      vis_compositing = 1;
      fluid_solver_and_vis_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  visRenderA (ColourPalette, &net, &sl);
	  
	  stability = lbmCycle (1, 1, 1, &lbm, &net);
	  
	  visRenderB (0, complete_image_name, ColourPalette, &net);
	  
	  // partial timings
	  elapsed_time = myClock () - fluid_solver_and_vis_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " FS + VIS, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      IsBenchSectionFinished (0.5, elapsed_time))
	    {
	      break;
	    }
	}
      fluid_solver_and_vis_time_steps = (int)(time_step * minutes / (3. * 0.5) - time_step);
      fluid_solver_and_vis_time = myClock ();
      
      for (time_step = 1; time_step <= fluid_solver_and_vis_time_steps; time_step++)
	{
	  visRenderA (ColourPalette, &net, &sl);
	  
	  stability = lbmCycle (1, 1, 1, &lbm, &net);
	  
	  visRenderB (0, complete_image_name, ColourPalette, &net);
	}
      fluid_solver_and_vis_time = myClock () - fluid_solver_and_vis_time;
      
      // benchmarking HemeLB's ray tracer without compositing
      
      vis_compositing = 0;
      vis_without_compositing_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  visRenderA (ColourPalette, &net, &sl);
	  
	  // partial timings
	  elapsed_time = myClock () - vis_without_compositing_time;
	  
	  if (time_step%bench_period == 1 && net.id == 0)
	    {
	      fprintf (stderr, " VIS - COMP, time: %.3f, time step: %i, time steps/s: %.3f\n",
		       elapsed_time, time_step, time_step / elapsed_time);
	    }
	  if (time_step%bench_period == 1 &&
	      IsBenchSectionFinished (0.5, elapsed_time))
	    {
	      break;
	    }
	}
      vis_without_compositing_time_steps = (int)(time_step * minutes / (3. * 0.5) - time_step);
      vis_without_compositing_time = myClock ();
      
      for (time_step = 1; time_step <= vis_without_compositing_time_steps; time_step++)
	{
	  visRenderA (ColourPalette, &net, &sl);
	}
      vis_without_compositing_time = myClock () - vis_without_compositing_time;
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
		   fluid_solver_time_steps / fluid_solver_time,
		   1.e-6 * lbm.total_fluid_sites / (fluid_solver_time / fluid_solver_time_steps),
		   fluid_solver_time);
	  
	  fprintf (timings_ptr, " FS + VIS, time steps per second: %.3f, time: %.3f\n\n",
		   fluid_solver_and_vis_time_steps / fluid_solver_and_vis_time, fluid_solver_and_vis_time);
	  
	  fprintf (timings_ptr, " VR - COMP, time steps per second: %.3f, time: %.3f\n\n",
		   vis_without_compositing_time_steps / vis_without_compositing_time, vis_without_compositing_time);
	}
    }
  
  if (net.id == 0)
    {
      fprintf (timings_ptr, "Opening output config file:\n %s\n\n", output_config_name);
      fflush (timings_ptr);
    }
  net.fo_time = myClock ();
  
  lbmWriteConfig (stability, output_config_name, &lbm, &net);
  
  net.fo_time = myClock () - net.fo_time;
  
  if (net.id == 0)
    {
      if (!is_bench)
	{
	  vis_pressure_min = lbmConvertPressureToPhysicalUnits (lbm_density_min * Cs2, &lbm);
	  vis_pressure_max = lbmConvertPressureToPhysicalUnits (lbm_density_max * Cs2, &lbm);
	  
	  vis_velocity_min = lbmConvertVelocityToPhysicalUnits (lbm_velocity_min, &lbm);
	  vis_velocity_max = lbmConvertVelocityToPhysicalUnits (lbm_velocity_max, &lbm);
	  
	  vis_stress_min = lbmConvertStressToPhysicalUnits (lbm_stress_min, &lbm);
	  vis_stress_max = lbmConvertStressToPhysicalUnits (lbm_stress_max, &lbm);
	  
	  fprintf (timings_ptr, "time steps per cycle: %i\n", lbm.period);
	  fprintf (timings_ptr, "pressure min, max (mmHg): %le, %le\n", vis_pressure_min, vis_pressure_max);
	  fprintf (timings_ptr, "velocity min, max (m/s) : %le, %le\n", vis_velocity_min, vis_velocity_max);
	  fprintf (timings_ptr, "stress   min, max (Pa)  : %le, %le\n", vis_stress_min, vis_stress_max);
	  fprintf (timings_ptr, "\n");
	  
	  for (int n = 0; n < lbm.inlets; n++)
	    {
	      fprintf (timings_ptr, "inlet id: %i, (flow rate)x(inlet area) (m^3/s): %le\n", n, lbm_inlet_flux[ n ]);
	    }
	  fprintf (timings_ptr, "\n");
	}
      fprintf (timings_ptr, "\n");
      fprintf (timings_ptr, "domain decomposition time (s):             %.3f\n", net.dd_time);
      fprintf (timings_ptr, "pre-processing buffer management time (s): %.3f\n", net.bm_time);
      fprintf (timings_ptr, "input configuration reading time (s):      %.3f\n", net.fr_time);
      fprintf (timings_ptr, "flow field outputting time (s):            %.3f\n", net.fo_time);
      
      total_time = myClock () - total_time;
      fprintf (timings_ptr, "total time (s):                            %.3f\n\n", total_time);
      
      if (net.procs > 1)
	{
	  fprintf (timings_ptr, "(vis+steering time) / fluid solver time: %.3e\n", steering_and_vis_time);
	}
      fprintf (timings_ptr, "Sub-domains info:\n\n");
      
      for (int n = 0; n < net.procs; n++)
	{
	  fprintf (timings_ptr, "rank: %i, fluid sites: %i\n", n, net.fluid_sites[ n ]);
	}
      
      fclose (timings_ptr);
    }
  
  visEnd (&sl);
  netEnd (&net);
  lbmEnd (&lbm);
  
  if (net.id == 0)
    {
      // there are some problems if the following function is called
      
      //pthread_join (network_thread, NULL);
      free(xdrSendBuffer_frame_details);
      free(xdrSendBuffer_pixel_data);
    }
  
  
#ifndef NOMPI
  net.err = MPI_Finalize ();
#endif
  
  return(0);
}

