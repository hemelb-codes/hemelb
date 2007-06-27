// In this file, the functions useful to initiate/end the LB simulation
// and perform the dynamics are reported

#include "config.h"


void lbmInit (char *system_file_name, char *parameters_file_name, char *checkpoint_file_name,
	      LBM *lbm, Net *net)
{
  // basically, this function call other ones only
  
  lbm->system_file_name     = system_file_name;
  lbm->parameters_file_name = parameters_file_name;
  lbm->checkpoint_file_name = checkpoint_file_name;
  
  lbmReadAndSetConfig (lbm, net);
  lbmReadParameters (lbm);
  lbmSetInitialConditions (lbm, net);
}


int lbmCycle (int write_checkpoint, int check_convergence, int *is_converged, LBM *lbm, Net *net)
{
  // the entire simulation time step takes place through this function
  
  double seconds;
  double f_eq[15];
  double omega;
  double density;
  double vx, vy, vz;
  double sum1, sum2;
  double stability_and_convergence_partial[3];
  double stability_and_convergence_total[3];
  double *f_old_p, *f_new_p;
  
  int i, l, m, n;
  int is_unstable;
  int *f_id_p;
  
  unsigned int site_data;
  
  short int *f_data_p;
  
  Velocity *vel_p;
  
  NeighProc *neigh_proc_p;
  
  
  is_unstable = 0;
  
  *is_converged = 0;
  
  sum1 = 0.0;
  sum2 = 0.0;
  
  omega = lbm->omega;
  
  seconds = myClock ();
  
  for (i = net->my_inner_sites;
       i < net->my_inner_sites + net->my_inter_sites; i++)
    {
      site_data = net->site_data[ i ];
      f_old_p = &f_old[ i*15 ];
      
      if (site_data == FLUID_TYPE)
	{
	  lbmFeq (f_old_p, &density, &vx, &vy, &vz, f_eq);
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_old_p[l] += omega * (f_old_p[l] - f_eq[l]);
	    }
	}
      else
	{
	  CalculateBC (f_old_p, site_data, &vx, &vy, &vz, lbm);
	}
      vel_p = &vel[ i ];
      
      sum1 += fabs(vel_p->x - vx) + fabs(vel_p->y - vy) + fabs(vel_p->z - vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel_p->x = vx;
      vel_p->y = vy;
      vel_p->z = vz;
    }
  net->timing[0] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  neigh_proc_p->f_to_send[ n ] = f_old[ neigh_proc_p->f_send_id[n] ];
	}
    }
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      neigh_proc_p = &net->inter_m_neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  neigh_proc_p->f_to_send[ n ] = f_old[ neigh_proc_p->f_send_id[n] ];
	}
    }
  net->timing[1] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      net->err = MPI_Isend (&neigh_proc_p->f_to_send[ 0 ],
			     neigh_proc_p->fs, MPI_DOUBLE,
			     neigh_proc_p->id, 10, MPI_COMM_WORLD,
			     &net->req[ 0 ][ net->id * net->procs + m ]);
      
      net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 0 ][ (net->id + net->procs) * net->procs + m ]);
    }
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      neigh_proc_p = &net->inter_m_neigh_proc[ m ];
      
      net->err = MPI_Issend (&neigh_proc_p->f_to_send[ 0 ],
			     neigh_proc_p->fs, MPI_DOUBLE,
			     neigh_proc_p->id, 10, MPI_COMM_WORLD,
			     &net->req[ 1 ][ net->id * net->procs + m ]);
      
      net->err = MPI_Irecv (&neigh_proc_p->f_to_recv[ 0 ],
			    neigh_proc_p->fs, MPI_DOUBLE,
			    neigh_proc_p->id, 10, MPI_COMM_WORLD,
			    &net->req[ 1 ][ (net->id + net->procs) * net->procs + m ]);
    }

  net->timing[2] += myClock () - seconds;
  seconds = myClock ();
  
  for (i = 0; i < net->my_inner_sites; i++)
    {
      site_data = net->site_data[ i ];
      f_old_p = &f_old[ i*15 ];
      
      if (site_data == FLUID_TYPE)
	{
	  lbmFeq (f_old_p, &density, &vx, &vy, &vz, f_eq);
	  
	  for (l = 0; l < 15; l++)
	    {
	      f_old_p[l] += omega * (f_old_p[l] - f_eq[l]);
	    }
	}
      else
	{
	  CalculateBC (f_old_p, site_data, &vx, &vy, &vz, lbm);
	}
      f_id_p = &f_id[ i*15 ];
      
      for (l = 0; l < 15; l++)
	{
	  if (f_old_p[l] < 0.) is_unstable = 1;
	  
	  f_new[ f_id_p[l] ] = f_old_p[l];
	}
      vel_p = &vel[ i ];
      
      sum1 += fabs(vel_p->x - vx) + fabs(vel_p->y - vy) + fabs(vel_p->z - vz);
      sum2 += fabs(vx) + fabs(vy) + fabs(vz);
      
      vel_p->x = vx;
      vel_p->y = vy;
      vel_p->z = vz;
    }
  net->timing[3] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      net->err = MPI_Wait (&net->req[ 1 ][ net->id * net->procs + m ], net->status);
      net->err = MPI_Wait (&net->req[ 1 ][ (net->id + net->procs) * net->procs + m ], net->status);
    }
  for (m = 0; m < net->neigh_procs; m++)
    {
      net->err = MPI_Wait (&net->req[ 0 ][ net->id * net->procs + m ], net->status);
      net->err = MPI_Wait (&net->req[ 0 ][ (net->id + net->procs) * net->procs + m ], net->status);
    }
  net->timing[2] += myClock () - seconds;
  seconds = myClock ();
  
  for (m = 0; m < net->inter_m_neigh_procs; m++)
    {
      neigh_proc_p = &net->inter_m_neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  f_new[ neigh_proc_p->f_recv_iv[n] ] = neigh_proc_p->f_to_recv[ n ];
	}
    }
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  f_new[ neigh_proc_p->f_recv_iv[n] ] = neigh_proc_p->f_to_recv[ n ];
	}
    }
  net->timing[1] += myClock () - seconds;
  seconds = myClock ();
  
  for (i = net->my_inner_sites;
       i < net->my_inner_sites + net->my_inter_sites; i++)
    {
      f_old_p = &f_old[ i*15 ];
      f_id_p = &f_id[ i*15 ];
      
      for (l = 0; l < 15; l++)
	{
	  if (f_old_p[l] < 0.) is_unstable = 1;
	  
	  f_new[ f_id_p[l] ] = f_old_p[l];
	}
    }
  f_old_p = f_old;
  f_old = f_new;
  f_new = f_old_p;
  
  net->timing[4] += myClock () - seconds;
  seconds = myClock ();
  
  free(lbm->inlet_density);
  free(lbm->outlet_density);
  
  lbmReadParameters (lbm);
  
  if (check_convergence)
    {
      if (net->procs > 1)
	{
	  stability_and_convergence_partial[ 0 ] = (double)is_unstable;
	  stability_and_convergence_partial[ 1 ] = sum1;
	  stability_and_convergence_partial[ 2 ] = sum2;
	  
	  net->err = MPI_Allreduce (stability_and_convergence_partial,
				    stability_and_convergence_total, 3,
				    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	  
	  sum1 = stability_and_convergence_total[ 1 ];
	  sum2 = stability_and_convergence_total[ 2 ];
	  
	  is_unstable = (stability_and_convergence_total[ 0 ] >= 1.);
	}
      
      if (sum1 <= sum2 * lbm->tolerance && sum2 > lbm->tolerance)
	{
	  *is_converged = 1;
	}
    }
  if (write_checkpoint)
    {
      lbmWriteConfig (!is_unstable, lbm->checkpoint_file_name, 1, lbm, net);
    }
  
  net->timing[5] += myClock () - seconds;
  
  if (is_unstable)
    {
      return UNSTABLE;
    }
  else
    {
      return STABLE;
    }
}

void lbmEnd (LBM *lbm)
{
  free(lbm->outlet_density);
  lbm->outlet_density = NULL;
  
  free(lbm->inlet_density);
  lbm->inlet_density = NULL;
  
  free(lbm->fluid_sites_per_block);
  lbm->fluid_sites_per_block = NULL;
}


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
