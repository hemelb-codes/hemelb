// In this file, the functions useful to initiate/end the LB simulation
// and perform the dynamics are reported

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
  
  int cycle_id;
  int total_time_steps, time_step, stability = STABLE;
  int depths;
  
  int fluid_solver_time_steps;
  int fluid_solver_and_vis_time_steps;
  int vis_without_compositing_time_steps;
  
  pthread_t network_thread;
  pthread_attr_t pthread_attrib;
  
  LBM lbm;
  
  Net net;
  
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
  else if (argc == 5)
    {
      is_bench = 0;
      lbm.cycles_max = atoi( argv[2] );
      lbm.period     = atoi( argv[3] );
      lbm.voxel_size = atof( argv[4] );
      
      if (lbm.cycles_max > 100)
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
  char output_image_name[256];
  char timings_name[256];
  char procs_string[256];
  char image_name[256];
  char output_directory[256];
  
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
  strcpy (output_directory, input_file_path);
  strcat (output_directory, "/Images/");
  mkdir  (output_directory, 0777);
  strcpy (output_image_name, output_directory);

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
  
  visInit (&net, &vis);
  
  visReadParameters (vis_parameters_name, &lbm, &net, &vis);
  
  
  if (!is_bench)
    {
      int is_finished = 0;
      
      total_time_steps = 0;
      
      simulation_time = myClock ();
      
      for (cycle_id = 0; cycle_id < lbm.cycles_max && !is_finished; cycle_id++)
	{
	  for (time_step = 0; time_step < lbm.period; time_step++)
	    {

//		usleep(40000);

	      int perform_rt = 0;
	      int write_image = 0;
	      int stream_image = 0; 
	      int lock_return = 0;
	      
	      total_time_steps++;
	      
	      if ((time_step + 1)%vis_image_freq == 0)
		{
		  // A write_image = 1;
		}

		  if (net.id == 0)
		    {
		      lock_return = pthread_mutex_trylock ( &LOCK );
		      printf("attempting to aquire mutex lock -> %i -> ", lock_return); 
		      if( lock_return == EBUSY ) { printf("lock busy\n"); } else { printf("aquired lock\n"); doRendering = 1; } printf("ShouldIRenderNow %i\n", ShouldIRenderNow); fflush(0x0);
		    }

	  UpdateSteerableParameters (&vis,&lbm);
	  MPI_Bcast (&doRendering, 1, MPI_INT, 0, MPI_COMM_WORLD);

	      if (total_time_steps%1 == 0)
		{

	 //    if(net.id==0) printf("cycle_id time_step setRendering doRendering %i %i %i %i\n", cycle_id, time_step, setRendering, doRendering); fflush(0x0); 

		  //UpdateSteerableParameters (&is_thread_locked, &vis);
		  
		 /* if (!is_thread_locked)
=======
		  if (is_thread_locked == 0)
>>>>>>> 1.31.2.9.2.28
		    {
<<<<<<< main.cc
		      stream_image = 1;
		    } */
		}

	      if (stream_image || write_image)
		{
		  doRendering = 1;
		}
	      
	      // Between the visRenderA/B calls, do not change any vis
	      // parameters.
	      
	//	if(setRendering==1) { doRendering = 1; setRendering = 0; }

	      if (doRendering)
		{
		  visRenderA (ColourPalette, &net);
		}

	      lbmVaryBoundaryDensities (cycle_id, time_step, &lbm);
	      
	      if (!check_conv)
		{
		  stability = lbmCycle (cycle_id, time_step, doRendering, &lbm, &net);
		}
	      else
		{
		  stability = lbmCycleConv (cycle_id, time_step, doRendering, &lbm, &net);
		}

	   /*   if (write_image)
		{
		  char time_step_string[256];
		  
		  // At this point output_image_name is appended with
		  // Images
		  strcpy ( image_name , output_image_name );

		  int time_steps = time_step + 1;
		  
		  while (time_steps < 100000000)
		    {
		      strcat ( image_name , "0" );
		      time_steps *= 10;
		    }
		  sprintf ( time_step_string, "%i", time_step + 1);
		  strcat ( image_name , time_step_string );
		  strcat ( image_name , ".dat" );
		} */

	      if (doRendering)
		{
		  visRenderB (write_image, image_name, ColourPalette, &net);
		}

	     /* if ((time_step + 1) % 10 == 0 ) {
                char snapshot_filename[255];
                snprintf(snapshot_filename, 255, "snapshot_%06i.bin", time_step+1);
                //if(net.id == 0) { printf("writing binary file %s....\n", snapshot_filename); fflush(NULL); }
                lbmWriteConfig (stability, snapshot_filename, &lbm, &net);
                //if(net.id == 0) { printf("done writing binary.\n"); fflush(NULL); }
              } */

	      if (net.id == 0)
		{
		if(doRendering==1) {
                  printf("sending signal to thread that frame is ready to go...\n"); fflush(0x0);
		  pthread_mutex_unlock (&LOCK);
		  pthread_cond_signal (&network_send_frame);
                  printf("...signal sent\n"); fflush(0x0);
                  doRendering=0;
		}
		} 
	      
	      if (stability == UNSTABLE)
		{
		  printf (" ATTENTION: INSTABILITY CONDITION OCCURRED\n");
		  printf (" AFTER %i total time steps\n", total_time_steps);
		  printf (" EXECUTION IS ABORTED\n");
#ifndef NOMPI
		  MPI_Abort (MPI_COMM_WORLD, 1);
#else
		  exit(1);
#endif
		  is_finished = 1;
		  break;
		}
	      if (stability == STABLE_AND_CONVERGED)
		{
		  is_finished = 1;
		  break;
		}
	      if (net.id == 0)
		{
		  //printf ("time step: %i\n", time_step+1);
		}
	    }
	  if (net.id == 0)
	    {
	      if (!check_conv)
		{
		  fprintf (timings_ptr, "cycle id: %i\n", cycle_id+1);
		  printf ("cycle id: %i\n", cycle_id+1);
		}
	      else
		{
		  fprintf (timings_ptr, "cycle id: %i, conv_error: %le\n", cycle_id+1, conv_error);
		  printf ("cycle id: %i, conv_error: %le\n", cycle_id+1, conv_error);
		}
                fflush(NULL);
	    }
	}
      simulation_time = myClock () - simulation_time;
      time_step = (1+min(time_step, lbm.period-1)) * min(cycle_id, lbm.cycles_max-1);
    }
  else // is_bench
    {
      double elapsed_time;
  
      int bench_period = (int)fmax(1., (1e+6 * net.procs) / lbm.total_fluid_sites);
      
      // benchmarking HemeLB's fluid solver only
      
      fluid_solver_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  stability = lbmCycle (0, 0, 0, &lbm, &net);
	  
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
	  stability = lbmCycle (0, 0, 1, &lbm, &net);
	}
      fluid_solver_time = myClock () - fluid_solver_time;
      
      
      // benchmarking HemeLB's fluid solver and ray tracer
      
      vis_image_freq = 1;
      vis_compositing = 1;
      fluid_solver_and_vis_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  visRenderA (ColourPalette, &net);
	  
	  stability = lbmCycle (0, 0, 1, &lbm, &net);
	  
	  visRenderB (0, image_name, ColourPalette, &net);
	  
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
	  visRenderA (ColourPalette, &net);
	  
	  stability = lbmCycle (0, 0, 1, &lbm, &net);
	  
	  visRenderB (0, image_name, ColourPalette, &net);
	}
      fluid_solver_and_vis_time = myClock () - fluid_solver_and_vis_time;
      
      // benchmarking HemeLB's ray tracer without compositing
      
      vis_compositing = 0;
      vis_without_compositing_time = myClock ();
      
      for (time_step = 1; time_step <= 1000000000; time_step++)
	{
	  visRenderA (ColourPalette, &net);
	  
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
	  visRenderA (ColourPalette, &net);
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
	  double pressure_min = lbmConvertPressureToPhysicalUnits (lbm_density_min * Cs2, &lbm);
	  double pressure_max = lbmConvertPressureToPhysicalUnits (lbm_density_max * Cs2, &lbm);
	  
	  double velocity_min = lbmConvertVelocityToPhysicalUnits (lbm_velocity_min, &lbm);
	  double velocity_max = lbmConvertVelocityToPhysicalUnits (lbm_velocity_max, &lbm);
	  
	  double stress_min = lbmConvertStressToPhysicalUnits (lbm_stress_min, &lbm);
	  double stress_max = lbmConvertStressToPhysicalUnits (lbm_stress_max, &lbm);
	  
	  fprintf (timings_ptr, "pressure min, max (mmHg): %le, %le\n", pressure_min, pressure_max);
	  fprintf (timings_ptr, "velocity min, max (m/s) : %le, %le\n", velocity_min, velocity_max);
	  fprintf (timings_ptr, "stress   min, max (Pa)  : %le, %le\n", stress_min, stress_max);
	}
      fprintf (timings_ptr, "\n");
      fprintf (timings_ptr, "domain decomposition time (s):             %.3f\n", net.dd_time);
      fprintf (timings_ptr, "pre-processing buffer management time (s): %.3f\n", net.bm_time);
      fprintf (timings_ptr, "input configuration reading time (s):      %.3f\n", net.fr_time);
      fprintf (timings_ptr, "flow field outputting time (s):            %.3f\n", net.fo_time);
      
      total_time = myClock () - total_time;
      fprintf (timings_ptr, "total time (s):                            %.3f\n\n", total_time);
      
      fprintf (timings_ptr, "Sub-domains info:\n\n");
      
      for (int n = 0; n < net.procs; n++)
	{
	  fprintf (timings_ptr, "rank: %i, fluid sites: %i\n", n, net.fluid_sites[ n ]);
	}
      
      fclose (timings_ptr);
    }
  
  visEnd ();
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

