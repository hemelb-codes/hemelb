// In this file the functions useful to discover the topology used and
// to create and delete the domain decomposition are reported
#include "config.h"


int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net)
{
  int i, j, k, ii, jj, kk;

  ProcBlock *proc_block_p;
  
  
  if (site_i < 0 || site_i >= sites_x ||
      site_j < 0 || site_j >= sites_y ||
      site_k < 0 || site_k >= sites_z)
    {
      return NULL;   // out of the bounding box
    }
  
  i = site_i >> shift;
  j = site_j >> shift;
  k = site_k >> shift;

  proc_block_p = &net->proc_block[(i * blocks_y + j) * blocks_z + k];

  if (proc_block_p->proc_id == NULL)
    {
      return NULL;   // an empty (solid) block is addressed
    }
  else
    {
      ii = site_i - (i << shift);
      jj = site_j - (j << shift);
      kk = site_k - (k << shift);
      
      return &proc_block_p->proc_id[(((ii << shift) + jj) << shift) + kk];
    }
}


unsigned int *netSiteMapPointer (int site_i, int site_j, int site_k, Net *net)
{
  // function useful to recover the site data through its specific
  // global coordinates (site_i,site_j,site_k)
  
  int i, j, k, ii, jj, kk;

  DataBlock *map_block_p;


  if (site_i < 0 || site_i >= sites_x ||
      site_j < 0 || site_j >= sites_y ||
      site_k < 0 || site_k >= sites_z)
    {
      return NULL;   // out of the bounding box
    }
  
  i = site_i >> shift;
  j = site_j >> shift;
  k = site_k >> shift;

  map_block_p = &net->map_block[(i * blocks_y + j) * blocks_z + k];

  if (map_block_p->site_data == NULL)
    {
      return NULL;   // an empty (solid) block is addressed
    }
  else
    {
      ii = site_i - (i << shift);
      jj = site_j - (j << shift);
      kk = site_k - (k << shift);
      
      return &map_block_p->site_data[(((ii << shift) + jj) << shift) + kk];
    }
}

//#undef MPICHX_TOPOLOGY_DEPTHS
#ifdef MPICHX_TOPOLOGY_DEPTHS

int netFindTopology (Net *net, int *depths)
{
  // the topology discovery mechanism is implemented in this
  // function.
  
  int *depth, **color;
  int machine_id, flag, is_found;
  int i, j, sum;
  
  *depths = 0;
  
  
  net->err = MPI_Attr_get (MPI_COMM_WORLD, MPICHX_TOPOLOGY_DEPTHS, &depth, &flag);
  
  if (net->err != MPI_SUCCESS || flag == 0)
    {
      return 0;
    }
  
  net->err = MPI_Attr_get (MPI_COMM_WORLD, MPICHX_TOPOLOGY_COLORS, &color, &flag);
  
  if (net->err != MPI_SUCCESS || flag == 0)
    {
      return 0;
    }

  net_machines = 0;
  
  net->machine_id = (int *)malloc(sizeof(int)* net->procs);
  net->procs_per_machine = (int *)malloc(sizeof(int)* net->procs);
  
  for (i = 0; i < net->procs; i++)
    {
      net->procs_per_machine[ i ] = 0;
    }
  for (i = 0; i < net->procs; i++)
    {
      if (depth[ i ] != 4) continue;

      *depths = max(*depths, depth[ i ]);
      
      for (j = 0, is_found = 0; j < net_machines && is_found == 0; j++)
	{
	  if (color[ i ][ 3 ] == net->machine_id[ j ])
	    {
	      is_found = 1;
	      ++net->procs_per_machine[ net->machine_id[j] ];
	    }
	}
      if (is_found == 1) continue;
      
      net->machine_id[ net_machines ] = color[ i ][ 3 ];
      ++net->procs_per_machine[ net_machines ];
      ++net_machines;
    }
  net_machines = max(1, net_machines);
  
  if (net_machines > MACHINES_MAX)
    {
      printf (" too many checked machines\n");
      printf (" the execution is terminated\n");
      net->err = MPI_Abort (MPI_COMM_WORLD, 1);
    }
  if (net_machines == 1)
    {
      for (i = 0; i < net->procs; i++)
	{
	  net->machine_id[ i ] = 0;
	}
        net->procs_per_machine[ 0 ] = net->procs;
    }
  else
    {
      for (i = 0; i < net->procs; i++)
	{
	  sum = 0;
	  machine_id = 0;
	  
	  is_found = 0;
	  
	  while (!is_found)
	    {
	      if (sum + net->procs_per_machine[ machine_id ] > i)
		{
		  is_found = 1;
		  continue;
		}
	      sum += net->procs_per_machine[ machine_id ];
	      ++machine_id;
	    }
	  net->machine_id[ i ] = machine_id;
	}
    }
  return 1;
}

#else

int netFindTopology (Net *net, int *depths)
{
  // the machine is assumed to be only one if this function is
  // used instead of the previous one

  *depths = 1;
  
  net_machines = 1;
  
  net->machine_id = (int *)malloc(sizeof(int) * net->procs);
  net->procs_per_machine = (int *)malloc(sizeof(int) * net_machines);
  
  for (int i = 0; i < net->procs; i++)
    {
      net->machine_id[ i ] = 0;
    }
  net->procs_per_machine[ 0 ] = net->procs;
  
  return 1;
}
#endif


void netInit (LBM *lbm, Net *net)
{
  // the domain partitioning technique and the management of the
  // buffers useful for the inter-processor communications are
  // implemented in this function
  
  double seconds;
  
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int i, j, k;
  int l, m, n;
  int mm;
  int sites_a, sites_b;
  int index_a;
  int sites_buffer_size;
  int unvisited_fluid_sites, partial_visited_fluid_sites;
  int proc_count;
  int fluid_sites_per_unit;
  int neigh_proc_index;
  int my_sites;
  int are_fluid_sites_incrementing;
  int is_inter_site, is_inner_site;
  int collision_offset[2][COLLISION_TYPES];
  int flag;
  int *proc_id_p;
  
  short int *f_data_p;
  
  unsigned int *site_data, *site_data_p;
  unsigned int site_map;
  
  bool *is_my_block;
  
  SiteLocation *site_location_a, *site_location_b;
  SiteLocation *site_location_a_p, *site_location_b_p;
  
  DataBlock *data_block_p;
  DataBlock *map_block_p;
  
  ProcBlock *proc_block_p;
  
  NeighProc *neigh_proc_p;
  
  
  net->fluid_sites = (int *)malloc(sizeof(int) * net->procs);
  
  sites_buffer_size = 10000;
  site_location_a = (SiteLocation *)malloc(sizeof(SiteLocation) * sites_buffer_size);
  site_location_b = (SiteLocation *)malloc(sizeof(SiteLocation) * sites_buffer_size);
  
  net->my_sites = 0;
  
  // a fast graph growing partitioning technique which spans the data
  // set only once is implemented here; the data set is explored by
  // means of the arrays "site_location_a[]" and "site_location_b[]"
  
  if (is_bench || net->procs == 1)
    {
      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites / (double)net->procs);
      proc_count = 0;
    }
  else
    {
      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites / (double)(net->procs - 1));
      proc_count = 1;
    }
  for (n = 0; n < net->procs; n++)
    {
      net->fluid_sites[ n ] = 0;
    }
  partial_visited_fluid_sites = 0;
  unvisited_fluid_sites = lbm->total_fluid_sites;
  
  seconds = myClock ();
  
  if (net_machines == 1 || net_machines == net->procs)
    {
      n = -1;
      
      for (i = 0; i < blocks_x; i++)
	{
	  for (j = 0; j < blocks_y; j++)
	    {
	      for (k = 0; k < blocks_z; k++)
		{
		  proc_block_p = &net->proc_block[ ++n ];
		  
		  if (proc_block_p->proc_id == NULL)
		    {
		      continue;
		    }
		  m = -1;
		  
		  for (site_i = i * block_size; site_i < i * block_size + block_size; site_i++)
		    {
		      for (site_j = j * block_size; site_j < j * block_size + block_size; site_j++)
			{
			  for (site_k = k * block_size; site_k < k * block_size + block_size; site_k++)
			    {
			      if (proc_block_p->proc_id[ ++m ] != -1)
				{
				  continue;
				}
			      proc_block_p->proc_id[ m ] = proc_count;
			      
			      if (proc_count == net->id)
				{
				  ++net->my_sites;
				}
			      ++partial_visited_fluid_sites;
			      
			      ++net->fluid_sites[ proc_count ];
			      
			      sites_a = 1;
			      site_location_a_p = &site_location_a[ 0 ];
			      site_location_a_p->i = site_i;
			      site_location_a_p->j = site_j;
			      site_location_a_p->k = site_k;
			      
			      are_fluid_sites_incrementing = 1;
			      
			      while (partial_visited_fluid_sites < fluid_sites_per_unit &&
				     are_fluid_sites_incrementing)
				{
				  sites_b = 0;
				  are_fluid_sites_incrementing = 0;
				  
				  for (index_a = 0;
				       index_a < sites_a && partial_visited_fluid_sites < fluid_sites_per_unit;
				       index_a++)
				    {
				      site_location_a_p = &site_location_a[ index_a ];
				      
				      for (l = 1;
					   l < 15 && partial_visited_fluid_sites < fluid_sites_per_unit;
					   l++)
					{
					  neigh_i = site_location_a_p->i + e_x[ l ];
					  neigh_j = site_location_a_p->j + e_y[ l ];
					  neigh_k = site_location_a_p->k + e_z[ l ];
					  
					  if (neigh_i == -1 || neigh_i == sites_x) continue;
					  if (neigh_j == -1 || neigh_j == sites_y) continue;
					  if (neigh_k == -1 || neigh_k == sites_z) continue;
					  
					  proc_id_p = netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
					  
					  if (proc_id_p == NULL || *proc_id_p != -1)
					    {
					      continue;
					    }
					  *proc_id_p = proc_count;
					  
					  ++partial_visited_fluid_sites;
					  
					  ++net->fluid_sites[ proc_count ];
					  
					  are_fluid_sites_incrementing = 1;
					  
					  if (sites_b == sites_buffer_size)
					    {
					      sites_buffer_size *= 2;
					      site_location_a = (SiteLocation *)realloc(site_location_a,
											sizeof(SiteLocation) * sites_buffer_size);
					      site_location_b = (SiteLocation *)realloc(site_location_b,
											sizeof(SiteLocation) * sites_buffer_size); 
					    }
					  site_location_b_p = &site_location_b[ sites_b ];
					  site_location_b_p->i = neigh_i;
					  site_location_b_p->j = neigh_j;
					  site_location_b_p->k = neigh_k;
					  ++sites_b;
					  
					  if (proc_count == net->id)
					    {
					      ++net->my_sites;
					    }
					}
				    }
				  site_location_a_p = site_location_a;
				  site_location_a = site_location_b;
				  site_location_b = site_location_a_p;
				  sites_a = sites_b;
				}
			      if (partial_visited_fluid_sites >= fluid_sites_per_unit)
				{
				  ++proc_count;
				  unvisited_fluid_sites -= partial_visited_fluid_sites;
				  fluid_sites_per_unit = (int)ceil((double)unvisited_fluid_sites / (double)(net->procs - proc_count));
				  partial_visited_fluid_sites = 0;
				}
			    }
			}
		    }
		}
	    }
	}
    }
  else
    {
      for (int unit_level = 1; unit_level >= 0; unit_level--)
	{
	  int up_units_max, marker, machine_id;
	  
	  if (unit_level == 1)
	    {
	      up_units_max = 1;
	      fluid_sites_per_unit = 0;
	    }
	  else
	    {
	      up_units_max = net_machines;
	      
	      if (is_bench)
		{
		  fluid_sites_per_unit = (int)ceil((double)unvisited_fluid_sites / (double)net->procs);
		  proc_count = 0;
		}
	      else
		{
		  fluid_sites_per_unit = (int)ceil((double)unvisited_fluid_sites / (double)(net->procs - 1));
		  proc_count = 1;
		}
	    }
	  for (int up_unit = 0; up_unit < up_units_max; up_unit++)
	    {
	      if (unit_level == 1)
		{
		  marker = -1;
		  
		  proc_count = net->procs;
		  
		  machine_id = proc_count - net->procs;
		  
		  if (is_bench)
		    {
		      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites *
						       (double)net->procs_per_machine[ machine_id ] / net->procs);
		    }
		  else
		    {
		      double weight;
		      
		      weight = (double)(net->procs_per_machine[machine_id] * net->procs) / (double)(net->procs - 1);
		      
		      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites * weight / net_machines);
		    }
		}
	      else
		{
		  marker = net->procs + up_unit;
		}
	      partial_visited_fluid_sites = 0;
	      
	      n = -1;
	      
	      for (i = 0; i < blocks_x; i++)
		{
		  for (j = 0; j < blocks_y; j++)
		    {
		      for (k = 0; k < blocks_z; k++)
			{
			  proc_block_p = &net->proc_block[ ++n ];
			  
			  if (proc_block_p->proc_id == NULL)
			    {
			      continue;
			    }
			  m = -1;
			  
			  for (site_i = i * block_size; site_i < i * block_size + block_size; site_i++)
			    {
			      for (site_j = j * block_size; site_j < j * block_size + block_size; site_j++)
				{
				  for (site_k = k * block_size; site_k < k * block_size + block_size; site_k++)
				    {
				      if (proc_block_p->proc_id[ ++m ] != marker)
					{
					  continue;
					}
				      proc_block_p->proc_id[ m ] = proc_count;
				      
				      if (proc_count == net->id)
					{
					  ++net->my_sites;
					}
				      ++partial_visited_fluid_sites;
				      
				      if (unit_level == 0) ++net->fluid_sites[ proc_count ];
				      
				      sites_a = 1;
				      site_location_a_p = &site_location_a[ 0 ];
				      site_location_a_p->i = site_i;
				      site_location_a_p->j = site_j;
				      site_location_a_p->k = site_k;
				      
				      are_fluid_sites_incrementing = 1;
				      
				      while (partial_visited_fluid_sites < fluid_sites_per_unit &&
					     are_fluid_sites_incrementing)
					{
					  sites_b = 0;
					  are_fluid_sites_incrementing = 0;
					  
					  for (index_a = 0;
					       index_a < sites_a && partial_visited_fluid_sites < fluid_sites_per_unit;
					       index_a++)
					    {
					      site_location_a_p = &site_location_a[ index_a ];
					      
					      for (l = 1;
						   l < 15 && partial_visited_fluid_sites < fluid_sites_per_unit;
						   l++)
						{
						  neigh_i = site_location_a_p->i + e_x[ l ];
						  neigh_j = site_location_a_p->j + e_y[ l ];
						  neigh_k = site_location_a_p->k + e_z[ l ];
						  
						  if (neigh_i == -1 || neigh_i == sites_x) continue;
						  if (neigh_j == -1 || neigh_j == sites_y) continue;
						  if (neigh_k == -1 || neigh_k == sites_z) continue;
						  
						  proc_id_p = netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
						  
						  if (proc_id_p == NULL || *proc_id_p != marker)
						    {
						      continue;
						    }
						  *proc_id_p = proc_count;
						  
						  ++partial_visited_fluid_sites;
						  
						  if (unit_level == 0) ++net->fluid_sites[ proc_count ];
						  
						  are_fluid_sites_incrementing = 1;
						  
						  if (sites_b == sites_buffer_size)
						    {
						      sites_buffer_size *= 2;
						      site_location_a = (SiteLocation *)realloc(site_location_a,
												sizeof(SiteLocation) * sites_buffer_size);
						      site_location_b = (SiteLocation *)realloc(site_location_b,
												sizeof(SiteLocation) * sites_buffer_size); 
						    }
						  site_location_b_p = &site_location_b[ sites_b ];
						  site_location_b_p->i = neigh_i;
						  site_location_b_p->j = neigh_j;
						  site_location_b_p->k = neigh_k;
						  ++sites_b;
						  
						  if (proc_count == net->id)
						    {
						      ++net->my_sites;
						    }
						}
					    }
					  site_location_a_p = site_location_a;
					  site_location_a = site_location_b;
					  site_location_b = site_location_a_p;
					  sites_a = sites_b;
					}
				      if (partial_visited_fluid_sites >= fluid_sites_per_unit)
					{
					  ++proc_count;
					  
					  if (unit_level == 0)
					    {
					      unvisited_fluid_sites -= partial_visited_fluid_sites;
					      fluid_sites_per_unit = (int)ceil((double)unvisited_fluid_sites / (double)(net->procs - proc_count));
					    }
					  partial_visited_fluid_sites = 0;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

  free(site_location_b);
  free(site_location_a);
  
  net->dd_time = myClock () - seconds;
  seconds = myClock ();
  
  // a map between the two-level data representation and the 1D
  // compact one is created here
  
  net->map_block = (DataBlock *)malloc(sizeof(DataBlock) * blocks);
  
  for (n = 0; n < blocks; n++)
    {
      net->map_block[ n ].site_data = NULL;
    }
  site_data = (unsigned int *)malloc(sizeof(unsigned int) * net->my_sites);
  
  is_my_block = (bool *)malloc(sizeof(bool) * blocks);
  
  for (n = 0; n < blocks; n++)
    {
      is_my_block[ n ] = 0;
    }
  my_sites = 0;
  
  for (n = 0; n < blocks; n++)
    {
      data_block_p = &net->data_block[ n ];
      
      if (data_block_p->site_data == NULL)
	{
	  continue;
	}
      proc_block_p = &net->proc_block[ n ];
      
      map_block_p = &net->map_block[ n ];
      map_block_p->site_data = (unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  if (proc_block_p->proc_id[ m ] == net->id)
	    {
	      if ((data_block_p->site_data[ m ] & SITE_TYPE_MASK) != SOLID_TYPE)
		{
		  map_block_p->site_data[ m ] = my_sites;
		  site_data[ my_sites ] = data_block_p->site_data[ m ];
		  ++my_sites;
		}
	      else
		{
		  map_block_p->site_data[ m ] = (1U << 31U);
		}
	      is_my_block[ n ] = 1;
	    }
	  else
	    {
	      map_block_p->site_data[ m ] = (1U << 31U);
	    }
	}
    }
  
  for (n = 0; n < blocks; n++)
    {
      if (net->data_block[ n ].site_data != NULL)
	{
	  free(net->data_block[ n ].site_data);
	  net->data_block[ n ].site_data = NULL;
	}
    }
  free(net->data_block);
  net->data_block = NULL;
  
  for (n = 0; n < blocks; n++)
    {
      if (!is_my_block[ n ])
	{
	  free(net->map_block[ n ].site_data);
	  net->map_block[ n ].site_data = NULL;
	}
    }
  free(is_my_block);
  
  // the number of inter- and intra-machine neighbouring processors,
  // interface-dependent and independent fluid sites and shared
  // distribution functions of the reference processor are calculated
  // here
  
  net->neigh_procs = 0;
  
  for (m = 0; m < NEIGHBOUR_PROCS_MAX; m++)
    {
      net->neigh_proc[ m ].fs = 0;
    }
  
  net->my_inter_sites = 0;
  net->my_inner_sites = 0;
  
  for (m = 0; m < COLLISION_TYPES; m++)
    {
      net->my_inter_collisions[ m ] = 0;
      net->my_inner_collisions[ m ] = 0;
    }
  net->shared_fs = 0;
  
  my_sites = 0;
  
  n = -1;
  
  for (i = 0; i < sites_x; i += block_size)
    {
      for (j = 0; j < sites_y; j += block_size)
	{
	  for (k = 0; k < sites_z; k += block_size)
	    {
	      map_block_p = &net->map_block[ ++n ];
	      
	      if (map_block_p->site_data == NULL) continue;
	      
	      proc_block_p = &net->proc_block[ n ];
	      
	      m = -1;
	      
	      for (site_i = i; site_i < i + block_size; site_i++)
		{
		  for (site_j = j; site_j < j + block_size; site_j++)
		    {
		      for (site_k = k; site_k < k + block_size; site_k++)
			{
			  if (proc_block_p->proc_id[ ++m ] != net->id)
			    {
			      continue;
			    }
			  is_inter_site = 0;
			  is_inner_site = 1;
			  
			  for (l = 1; l < 15; l++)
			    {
			      neigh_i = site_i + e_x[ l ];
			      neigh_j = site_j + e_y[ l ];
			      neigh_k = site_k + e_z[ l ];
			      
			      proc_id_p = netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
			      
			      if (proc_id_p == NULL || *proc_id_p == net->id || *proc_id_p == (1 << 30))
				{
				  continue;
				}
			      is_inner_site = 0;
			      is_inter_site = 1;
			      
			      for (mm = 0, flag = 1; mm < net->neigh_procs && flag; mm++)
				{
				  neigh_proc_p = &net->neigh_proc[ mm ];
				  
				  if (*proc_id_p == neigh_proc_p->id)
				    {
				      flag = 0;
				      ++neigh_proc_p->fs;
				      ++net->shared_fs;
				    }
				}
			      if (flag)
				{
				  if (net->neigh_procs == NEIGHBOUR_PROCS_MAX)
				    {
				      printf (" too many intra machine, inter processor neighbours\n");
				      printf (" the execution is terminated\n");
#ifndef NOMPI
				      net->err = MPI_Abort (MPI_COMM_WORLD, 1);
#else
				      exit(1);
#endif
				    }
				  neigh_proc_p = &net->neigh_proc[ net->neigh_procs ];
				  neigh_proc_p->id = *proc_id_p;
				  ++neigh_proc_p->fs;
				  ++net->neigh_procs;
				  ++net->shared_fs;
				}
			    }
			  l = lbmCollisionType (site_data[ my_sites ]);
			  ++my_sites;
			  
			  if (is_inner_site)
			    {
			      ++net->my_inner_sites;
			      
			      if (l == 0)
				{
				  map_block_p->site_data[ m ] = net->my_inner_collisions[l];
				}
			      else
				{
				  map_block_p->site_data[ m ] = 50000000 * (10 + (l-1)) + net->my_inner_collisions[l];
				}
			      ++net->my_inner_collisions[ l ];
			    }
			  else if (is_inter_site)
			    {
			      ++net->my_inter_sites;
			      
			      if (l == 0)
				{
				  map_block_p->site_data[ m ] = 1000000000 + net->my_inter_collisions[l];
				}
			      else
				{
				  map_block_p->site_data[ m ] = 50000000 * (20 + l) + net->my_inter_collisions[l];
				}
			      ++net->my_inter_collisions[ l ];
			    }
			}
		    }
		}
	    }
	}
    }
  
  collision_offset[0][0] = 0;
  
  for (l = 1; l < COLLISION_TYPES; l++)
    {
      collision_offset[0][l] = collision_offset[0][l-1] + net->my_inner_collisions[l-1];
    }
  collision_offset[1][0] = net->my_inner_sites;
  
  for (l = 1; l < COLLISION_TYPES; l++)
    {
      collision_offset[1][l] = collision_offset[1][l-1] + net->my_inter_collisions[l-1];
    }
  for (n = 0; n < blocks; n++)
    {
      map_block_p = &net->map_block[ n ];
      
      if (map_block_p->site_data == NULL) continue;
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  site_data_p = &map_block_p->site_data[ m ];
	  
	  if (*site_data_p & (1U << 31U)) continue;
	  
	  if (*site_data_p < 500000000)
	    {
	      continue;
	    }
	  for (l = 1; l < COLLISION_TYPES; l++)
	    {
	      if (*site_data_p >= 50000000 * (10 + (l-1)) &&
		  *site_data_p  < 50000000 * (10 +  l))
		{
		  *site_data_p += collision_offset[0][l] - 50000000 * (10 + (l-1));
		  break;
		}
	    }
	  for (l = 0; l < COLLISION_TYPES; l++)
	    {
	      if (*site_data_p >= 50000000 * (20 +  l) &&
		  *site_data_p  < 50000000 * (20 + (l+1)))
		{
		  *site_data_p += collision_offset[1][l] - 50000000 * (20 + l);
		  break;
		}
	    }
	}
    }
  
  f_old = (double *)malloc(sizeof(double) * (net->my_sites * 15 + 1 + net->shared_fs));
  f_new = (double *)malloc(sizeof(double) * (net->my_sites * 15 + 1 + net->shared_fs));
  
  // the precise interface-dependent data (interface-dependent fluid
  // site locations and identifiers of the distribution functions
  // streamed between different partitions) are collected and the
  // buffers needed for the communications are set from here
  
  f_data = (short int *)malloc(sizeof(short int) * 4 * net->shared_fs);
  
  f_recv_iv = (int *)malloc(sizeof(int) * net->shared_fs);
  
  net->shared_fs = 0;
  
  for (n = 0; n < net->neigh_procs; n++)
    {
      net->neigh_proc[ n ].f_data    = &f_data[ net->shared_fs<<2 ];
      net->neigh_proc[ n ].f_head    = net->my_sites * 15 + 1 + net->shared_fs;
      net->neigh_proc[ n ].f_recv_iv = &f_recv_iv[ net->shared_fs ];
      
      net->shared_fs += net->neigh_proc[ n ].fs;
      net->neigh_proc[ n ].fs = 0;
    }
  
  if (net->my_sites > 0)
    {
      f_id = (int *)malloc(sizeof(int) * (net->my_sites * 15));
      
      net_site_data = (unsigned int *)malloc(sizeof(unsigned int) * net->my_sites);
    }
  
  net->from_proc_id_to_neigh_proc_index = (short int *)malloc(sizeof(short int) * net->procs);
  
  for (m = 0; m < net->procs; m++)
    {
      net->from_proc_id_to_neigh_proc_index[ m ] = -1;
    }
  for (m = 0; m < net->neigh_procs; m++)
    {
      net->from_proc_id_to_neigh_proc_index[ net->neigh_proc[m].id ] = m;
    }
  
  my_sites = 0;
  
  n = -1;
  
  for (i = 0; i < sites_x; i += block_size)
    {
      for (j = 0; j < sites_y; j += block_size)
	{
	  for (k = 0; k < sites_z; k += block_size)
	    {
	      map_block_p = &net->map_block[ ++n ];
	      
	      if (map_block_p->site_data == NULL) continue;

	      proc_block_p = &net->proc_block[ n ];
	      
	      m = -1;
	      
	      for (site_i = i; site_i < i + block_size; site_i++)
		{
		  for (site_j = j; site_j < j + block_size; site_j++)
		    {
		      for (site_k = k; site_k < k + block_size; site_k++)
			{
			  if (proc_block_p->proc_id[ ++m ] != net->id)
			    {
			      continue;
			    }
			  site_map = map_block_p->site_data[ m ];
			  
			  f_id[ site_map*15+0 ] = site_map * 15 + 0;
			  
			  for (l = 1; l < 15; l++)
			    {
			      neigh_i = site_i + e_x[ l ];
			      neigh_j = site_j + e_y[ l ];
			      neigh_k = site_k + e_z[ l ];
			      
			      f_id[ site_map*15+l   ] = net->my_sites * 15;
			      
			      proc_id_p = netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
			      
			      if (proc_id_p == NULL || *proc_id_p == 1 << 30)
				{
				  continue;
				}
			      site_data_p = netSiteMapPointer (neigh_i, neigh_j, neigh_k, net);
			      
			      if (*proc_id_p == net->id)
				{
				  f_id[ site_map*15+l   ] = *site_data_p * 15 + l;
				  continue;
				}
			      neigh_proc_index = net->from_proc_id_to_neigh_proc_index[ *proc_id_p ];
			      
			      neigh_proc_p = &net->neigh_proc[ neigh_proc_index ];
			      
			      f_data_p = &neigh_proc_p->f_data[ neigh_proc_p->fs<<2 ];
			      f_data_p[ 0 ] = site_i;
			      f_data_p[ 1 ] = site_j;
			      f_data_p[ 2 ] = site_k;
			      f_data_p[ 3 ] = l;
			      ++neigh_proc_p->fs;
			    }
			  net_site_data[ site_map ] = site_data[ my_sites ];
			  ++my_sites;
			}
		    }
		}
	    }
	}
    }
  free(site_data);
  
  // point-to-point communications are performed to match data to be
  // sent to/receive from different partitions; in this way, the
  // communication of the locations of the interface-dependent fluid
  // sites and the identifiers of the distribution functions which
  // propagate to different partitions is avoided (only their values
  // will be communicated)
  
#ifndef NOMPI
  net->req = (MPI_Request **)malloc(sizeof(MPI_Request *) * COMMS_LEVELS);
  
  for (m = 0; m < COMMS_LEVELS; m++)
    {
      net->req[ m ] = (MPI_Request *)malloc(sizeof(MPI_Request) * (2 * net->procs));
    }
#endif
  
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      if (neigh_proc_p->id > net->id)
  	{
#ifndef NOMPI
	  net->err = MPI_Isend (&neigh_proc_p->f_data[ 0 ],
				neigh_proc_p->fs * 4, MPI_SHORT,
				neigh_proc_p->id, 10, MPI_COMM_WORLD,
				&net->req[ 0 ][ m ]);
#endif
  	}
      else
  	{
#ifndef NOMPI
  	  net->err = MPI_Irecv (&neigh_proc_p->f_data[ 0 ],
  				neigh_proc_p->fs * 4, MPI_SHORT,
  				neigh_proc_p->id, 10, MPI_COMM_WORLD,
  				&net->req[ 0 ][ net->neigh_procs + m ]);
#endif
  	}
    }
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      if (neigh_proc_p->id > net->id)
  	{
#ifndef NOMPI
  	  net->err = MPI_Wait (&net->req[ 0 ][ m ], net->status);
#endif
  	}
      else
  	{
#ifndef NOMPI
  	  net->err = MPI_Wait (&net->req[ 0 ][ net->neigh_procs + m ], net->status);
#endif
  	  for (n = 0; n < neigh_proc_p->fs*4; n += 4)
  	    {
	      f_data_p = &neigh_proc_p->f_data[ n ];
	      
  	      l = f_data_p[ 3 ];
  	      f_data_p[ 0 ] += e_x[ l ];
  	      f_data_p[ 1 ] += e_y[ l ];
  	      f_data_p[ 2 ] += e_z[ l ];
  	      f_data_p[ 3 ] = inv_dir[ l ];
  	    }
  	}
    }
  
  int f_count = net->my_sites * 15;
 
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  f_data_p = &neigh_proc_p->f_data[ n*4 ];
	  i = f_data_p[ 0 ];
	  j = f_data_p[ 1 ];
	  k = f_data_p[ 2 ];
	  l = f_data_p[ 3 ];
	  
	  site_map = *netSiteMapPointer (i, j, k, net);
	  
	  f_id[ site_map * 15 + l ] = ++f_count;
	  neigh_proc_p->f_recv_iv[ n ] = site_map * 15 + inv_dir[ l ];
	}
    }
  free(f_data);
  
  net->bm_time = myClock () - seconds;
}


void netEnd (Net *net)
{
  // the allocated data are freed with this function
  
  int i;
  
  
  free(net->from_proc_id_to_neigh_proc_index);
  net->from_proc_id_to_neigh_proc_index = NULL;
  
  for (i = 0; i < blocks; i++)
    {
      if (net->map_block[ i ].site_data != NULL)
	{
	  free(net->map_block[ i ].site_data);
	  net->map_block[ i ].site_data = NULL;
	}
    }
  
  free(f_recv_iv);
  
  
  free(net->map_block);
  net->map_block = NULL;
  
  for (i = 0; i < blocks; i++)
    {
      if (net->proc_block[ i ].proc_id != NULL)
	{
	  free(net->proc_block[ i ].proc_id);
	  net->proc_block[ i ].proc_id = NULL;
	}
    }
  free(net->proc_block);
  net->proc_block = NULL;
  
  if (net->my_sites > 0)
    {
      free(net_site_data);
      net_site_data = NULL;
      
      free(f_id);
      f_id = NULL;
    }
  free(f_new);
  f_new = NULL;
  free(f_old);
  f_old = NULL;
  
#ifndef NOMPI
  for (i = 0; i < COMMS_LEVELS; i++)
    {
      free(net->req[ i ]);
      net->req[ i ] = NULL;
    }
  free(net->req);
  net->req = NULL;
#endif
  
  free(net->fluid_sites);
  
  free(net->procs_per_machine);
  free(net->machine_id);
}

