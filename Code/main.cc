// In this file, the functions useful to initiate/end the LB simulation
// and perform the dynamics are reported

#include "config.h"


void usage (char *progname)
{
  printf ("Usage: %s input file name with config, parameters, output, checkpoint file names\n", progname);
}


int main (int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output
  
  double seconds;
  double total_time;
  
  int time_step, stability, is_converged;
  int write_checkpoint, do_convergence;
  int n;
  int checkpoint_count = 0;
  int convergence_count = 0;
  int required_args = 2;

  LBM lbm;
  
  Net net;
  
  
  net.err = MPI_Init (&argc, &argv);
  net.err = MPI_Comm_size (MPI_COMM_WORLD, &net.procs);
  net.err = MPI_Comm_rank (MPI_COMM_WORLD, &net.id);

  if (argc != required_args)
    {
      if (net.id == 0) usage(argv[0]);
      
      net.err = MPI_Abort (MPI_COMM_WORLD, 1);
      net.err = MPI_Finalize ();
    }
  
  char *input_file_name(argv[1]);
  char input_config_name[80];
  char input_parameters_name[80];
  char output_config_name[80];
  char checkpoint_config_name[80];
  
  FILE *input_file = fopen (input_file_name, "r");
  
  
  fscanf (input_file, "%s ", input_config_name);
  fscanf (input_file, "%s ", input_parameters_name);
  fscanf (input_file, "%s ", output_config_name);
  fscanf (input_file, "%s ", checkpoint_config_name);
  
  fclose (input_file);
  
  if (net.id == 0)
    {
      printf ("\n");
      printf ("***********************************************************\n");
      printf("Opening config file: %s\n", input_config_name);
      printf("Opening parameters file: %s\n", input_parameters_name);
    }
  
  lbmInit (input_config_name, input_parameters_name, checkpoint_config_name, &lbm, &net);
  
  for (n = 0; n < 7; n++) net.timing[ n ] = 0.;
  
  seconds = myClock ();
  
  checkpoint_count = 0;
  convergence_count = 0;
  
  for (time_step = 1; time_step <= lbm.time_steps_max; time_step++)
    {
      write_checkpoint = 0;
      do_convergence = 0;
      
      if (++checkpoint_count >= lbm.checkpoint_frequency)
	{
	  write_checkpoint = 1;
	  checkpoint_count = 0;
	}
      if (++convergence_count >= lbm.convergence_frequency)
	{
	  do_convergence = 1;
	  convergence_count = 0;
	}
      stability = lbmCycle (write_checkpoint, do_convergence, &is_converged, &lbm, &net);
      
      if (stability == UNSTABLE || is_converged) break;
    }
  time_step = min(time_step, lbm.time_steps_max);
  
  seconds = myClock () - seconds;
  
  if (net.id == 0)
    {
      if (stability == STABLE)
  	{
  	  if (!is_converged)
  	    {
  	      printf (" ATTENTION: SIMULATION NOT CONVERGED\n");
  	    }
  	  
  	  printf (" fluid sites: %4i, MLSUPS: %.3f, time steps: %i \n procs: %i, machines: %i\n",
  		  lbm.total_fluid_sites,
  		  time_step * 1.E-6 * lbm.total_fluid_sites / seconds, time_step,
		  net.procs, net.machines);
  	}
      else
  	{
  	  printf (" ATTENTION: INSTABILITY CONDITION OCCURRED\n");
  	  printf (" AFTER %i time steps\n", time_step);
  	}
    }
  
  seconds = myClock ();
  
  lbmWriteConfig (stability, output_config_name, 0, &lbm, &net);
  
  net.timing[6] = myClock () - seconds;
  
  net.err = MPI_Barrier (MPI_COMM_WORLD);
  
  if (net.id == 0)
    {
      printf ("timings results (seconds/time step)\n");
      printf ("1) rank\n");
      printf ("2) collisions at inter-machine lattice sites\n");
      printf ("3) copying to the buffers to send/recv\n");
      printf ("4) communicational time\n");
      printf ("5) collisions + streaming from inner sites\n");
      printf ("6) streaming from inter-machine lattice sites\n");
      printf ("7) all-reduce + convergence test + checkpoint\n");
      printf ("8) output (1 time only, seconds) \n");
    }
  total_time = 0.;
  
  for (n = 0; n < 6; n++) total_time += net.timing[ n ];
  
  for (n = 0; n < 6; n++)
    {
      if (net.timing[ n ] < 1.e-5 * total_time)
	{
	  net.timing[ n ] = 0.;
	}
      else
	{
	  net.timing[ n ] /= time_step;
	}
    }
  
  printf ("%i %.3e %.3e %.3e %.3e %.3e %.3e %.3e 1) -> 8)\n", net.id,
	  net.timing[0], net.timing[1], net.timing[2],
	  net.timing[3], net.timing[4], net.timing[5], net.timing[6]);
  
  lbmEnd (&lbm);
  netEnd (&net);
  
  net.err = MPI_Finalize ();
  
  return(0);
}
