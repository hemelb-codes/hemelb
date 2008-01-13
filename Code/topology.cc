// In this file the functions useful to discover the topology used and
// to create and delete the domain decomposition are reported
#include "config.h"


short int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net)
{
  // function useful to recover the processor rank corresponding to a
  // specific lattice site (with global coordinates
  // (site_i,site_j,site_k)
  
  if (site_i < 0 || site_i >= sites_x ||
      site_j < 0 || site_j >= sites_y ||
      site_k < 0 || site_k >= sites_z)
    {
      return NULL;
    }

  int i = site_i >> shift;
  int j = site_j >> shift;
  int k = site_k >> shift;
  
  return &net->proc_id[ (i * blocks_y + j) * blocks_z + k ];
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

///#undef MPICHX_TOPOLOGY_DEPTHS
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

  net->machines = 0;
  
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
      
      for (j = 0, is_found = 0; j < net->machines && is_found == 0; j++)
	{
	  if (color[ i ][ 3 ] == net->machine_id[ j ])
	    {
	      is_found = 1;
	      ++net->procs_per_machine[ net->machine_id[j] ];
	    }
	}
      if (is_found == 1) continue;
      
      net->machine_id[ net->machines ] = color[ i ][ 3 ];
      ++net->procs_per_machine[ net->machines ];
      ++net->machines;
    }
  net->machines = max(1, net->machines);
  
  if (net->machines > MACHINES_MAX)
    {
      printf (" too many checked machines\n");
      printf (" the execution is terminated\n");
      net->err = MPI_Abort (MPI_COMM_WORLD, 1);
    }
  if (net->machines == 1)
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
  
  net->machines = 1;
  
  net->machine_id = (int *)malloc(sizeof(int) * net->procs);
  net->procs_per_machine = (int *)malloc(sizeof(int) * net->machines);
  
  for (int i = 0; i < net->procs; i++)
    {
      net->machine_id[ i ] = 0;
    }
  net->procs_per_machine[ 0 ] = net->procs;
  
  return 1;
}

#endif


void netInit (LBM *lbm, Net *net, Vis *vis, int proc_sites[])
{
  // the domain partitioning technique and the management of the
  // buffers useful for the inter-processor communications are
  // implemented in this function
  
  XDR xdr_config;
  
  double seconds;
  
  int n_x[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, +0, +0, +0, +0, +0, +0, +0, +0, +1, +1, +1, +1, +1, +1, +1, +1, +1};
  int n_y[] = {-1, -1, -1, +0, +0, +0, +1, +1, +1, -1, -1, -1, +0, +0, +1, +1, +1, -1, -1, -1, +0, +0, +0, +1, +1, +1};  
  int n_z[] = {-1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1};
  int site_i, site_j, site_k;
  int neigh_i, neigh_j, neigh_k;
  int block_id;
  int i, j, k;
  int l, m, n;
  int mm;
  int blocks_a, blocks_b;
  int index_a;
  int blocks_buffer_size;
  int partial_visited_fluid_sites;
  int proc_count;
  int fluid_sites_per_unit;
  int unit_level, up_units_max, up_unit, marker;
  int machine_id;
  int neigh_proc_index;
  int my_sites;
  int are_fluid_sites_incrementing;
  int is_inter_site, is_inner_site;
  int collision_offset[2][COLLISION_TYPES];
  int flag;
  int use_fast_dd;
  int clusters_max;
  int dummy = 0;
  
  short int proc_id, *proc_id_p, neigh_proc_id;
  short int *f_data_p;

  unsigned int *site_data, *site_data_p;
  unsigned int site_map;
  
  unsigned short int cluster_block_max_i, cluster_block_max_j, cluster_block_max_k;
  
  BlockLocation *block_location_a, *block_location_b;
  BlockLocation *block_location_a_p, *block_location_b_p;
  
  DataBlock *data_block_p;
  DataBlock *map_block_p;
  
  NeighProc *neigh_proc_p;
  
  Cluster *cluster_p;
  
  
  
  if (net->machines == 1 || net->machines == net->procs)
    {
      // a (fast) single-level domain decomposition method will be employed
      use_fast_dd = 1;
    }
  else
    {
      // a (slower) two-lelvel version of the domain decomposition method will be employed
      use_fast_dd = 0;
    }
  
  seconds = myClock ();
  
  net->proc_id = (short int *)malloc(sizeof(short int) * blocks);
  
  clusters_max = 20;
  vis->clusters = 0;
  vis->cluster = (Cluster *)malloc(sizeof(Cluster) * clusters_max);
  
  marker = -1;
  
  for (n = 0; n < blocks; n++)
    {
      if (lbm->fluid_sites_per_block[ n ] == 0)
	{
	  net->proc_id[ n ] = 1 << 14;
	}
      else
	{
	  net->proc_id[ n ] = marker;
	}
    }
  
  blocks_buffer_size = 10000;
  
  block_location_a = (BlockLocation *)malloc(sizeof(BlockLocation) * blocks_buffer_size);
  block_location_b = (BlockLocation *)malloc(sizeof(BlockLocation) * blocks_buffer_size);
  
  cluster_p = NULL;
  cluster_block_max_i = dummy;
  cluster_block_max_j = dummy;
  cluster_block_max_k = dummy;
  
  if (use_fast_dd)
    {
      // a fast graph growing partitioning technique which spans the
      // data set only once is implemented here; the data set is
      // explored at the coarse level of the data hierarchy by means
      // of the arrays "block_location_a[]" and "block_location_b[]"
      
#ifndef FREEROOT
      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites / (double)net->procs);
      proc_count = 0;
#else
      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites / (double)(net->procs - 1));
      proc_count = 1;
#endif
      
      partial_visited_fluid_sites = 0;
      
      n = -1;
      
      for (i = 0; i < blocks_x; i++)
	{
	  for (j = 0; j < blocks_y; j++)
	    {
	      for (k = 0; k < blocks_z; k++)
		{
		  if (*(proc_id_p = &net->proc_id[ ++n ]) != marker) continue;
		  
		  *proc_id_p = proc_count;
		  
		  if (proc_count == net->id)
		    {
		      if (vis->clusters == clusters_max)
			{
			  clusters_max <<= 1;
			  vis->cluster = (Cluster *)realloc(vis->cluster, sizeof(Cluster) * clusters_max);
			}
		      cluster_p = &vis->cluster[ vis->clusters ];
		      ++vis->clusters;
		      
		      cluster_p->block_min[0] = i;
		      cluster_p->block_min[1] = j;
		      cluster_p->block_min[2] = k;
		      
		      cluster_block_max_i = i;
		      cluster_block_max_j = j;
		      cluster_block_max_k = k;
		    }
		  partial_visited_fluid_sites += lbm->fluid_sites_per_block[ n ];
		  
		  blocks_a = 1;
		  block_location_a_p = &block_location_a[ 0 ];
		  block_location_a_p->i = i;
		  block_location_a_p->j = j;
		  block_location_a_p->k = k;
		  
		  are_fluid_sites_incrementing = 1;
		  
		  while (partial_visited_fluid_sites < fluid_sites_per_unit &&
			 are_fluid_sites_incrementing)
		    {
		      blocks_b = 0;
		      are_fluid_sites_incrementing = 0;
		      
		      for (index_a = 0;
			   index_a < blocks_a && partial_visited_fluid_sites < fluid_sites_per_unit;
			   index_a++)
			{
			  block_location_a_p = &block_location_a[ index_a ];
			  
			  for (l = 0;
			       l < 26 && partial_visited_fluid_sites < fluid_sites_per_unit;
			       l++)
			    {
			      neigh_i = block_location_a_p->i + n_x[ l ];
			      neigh_j = block_location_a_p->j + n_y[ l ];
			      neigh_k = block_location_a_p->k + n_z[ l ];
			      
			      if (neigh_i == -1 || neigh_i == blocks_x) continue;
			      if (neigh_j == -1 || neigh_j == blocks_y) continue;
			      if (neigh_k == -1 || neigh_k == blocks_z) continue;
			      
			      block_id = (neigh_i * blocks_y + neigh_j) * blocks_z + neigh_k;
			      
			      if (*(proc_id_p = &net->proc_id[ block_id ]) != marker)
				{
				  continue;
				}
			      *proc_id_p = proc_count;
			      
			      partial_visited_fluid_sites += lbm->fluid_sites_per_block[ block_id ];
			      
			      are_fluid_sites_incrementing = 1;
			      
			      if (blocks_b == blocks_buffer_size)
				{
				  blocks_buffer_size *= 2;
				  block_location_a = (BlockLocation *)realloc(block_location_a,
									      sizeof(BlockLocation) * blocks_buffer_size);
				  block_location_b = (BlockLocation *)realloc(block_location_b,
									      sizeof(BlockLocation) * blocks_buffer_size); 
				}
			      block_location_b_p = &block_location_b[ blocks_b ];
			      block_location_b_p->i = neigh_i;
			      block_location_b_p->j = neigh_j;
			      block_location_b_p->k = neigh_k;
			      ++blocks_b;
			      
			      if (proc_count == net->id)
				{
				  cluster_p->block_min[0] = min((int)cluster_p->block_min[0], neigh_i);
				  cluster_p->block_min[1] = min((int)cluster_p->block_min[1], neigh_j);
				  cluster_p->block_min[2] = min((int)cluster_p->block_min[2], neigh_k);
				  
				  cluster_block_max_i = max((int)cluster_block_max_i, neigh_i);
				  cluster_block_max_j = max((int)cluster_block_max_j, neigh_j);
				  cluster_block_max_k = max((int)cluster_block_max_k, neigh_k);
				}
			    }
			}
		      
		      if (are_fluid_sites_incrementing)
			{
			  block_location_a_p = block_location_a;
			  block_location_a = block_location_b;
			  block_location_b = block_location_a_p;
			  blocks_a = blocks_b;
			}
		    }
		  if (net->id == proc_count)
		    {
		      cluster_p->x[0] = cluster_p->block_min[0] * block_size - 0.5F * sites_x;
		      cluster_p->x[1] = cluster_p->block_min[1] * block_size - 0.5F * sites_y;
		      cluster_p->x[2] = cluster_p->block_min[2] * block_size - 0.5F * sites_z;
		      
		      cluster_p->blocks_x = 1 + cluster_block_max_i - cluster_p->block_min[0];
		      cluster_p->blocks_y = 1 + cluster_block_max_j - cluster_p->block_min[1];
		      cluster_p->blocks_z = 1 + cluster_block_max_k - cluster_p->block_min[2];
		    }
		  if (partial_visited_fluid_sites >= fluid_sites_per_unit)
		    {
		      ++proc_count;
		      partial_visited_fluid_sites = 0;
		    }
		}
	    }
	}
    }
  else
    {
      // a slower graph growing partitioning technique which spans the
      // data set 2xM times where M is the number of machines used is
      // implemented here; the data set is explored at the coarse
      // level of the data hierarchy by means of the arrays
      // "block_location_a[]" and "block_location_b[]"
      
      for (unit_level = 1; unit_level >= 0; unit_level--)
	{
	  if (unit_level == 1)
	    {
	      up_units_max = 1;
	      fluid_sites_per_unit = dummy;
	    }
	  else
	    {
	      up_units_max = net->machines;
	      
//#ifndef FREEROOT
	      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites / (double)net->procs);
//#else
//	      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites / (double)(net->procs - 1));
//#endif
	    }
	  
	  for (up_unit = 0; up_unit < up_units_max; up_unit++)
	    {
	      if (unit_level == 1)
		{
		  marker = -1;
		  
//#ifndef FREEROOT
		  proc_count = net->procs;
//#else
//		  proc_count = net->procs - 1;
//#endif
		}
	      else
		{
//#ifndef FREEROOT
		  marker = net->procs + up_unit;
//#else
//		  marker = net->procs - 1 + up_unit;
//#endif
		  
		  proc_count = 0;

		  for (n = 0; n < up_unit; n++)
		    {
		      proc_count += net->procs_per_machine[ n ];
		    }
//#ifdef FREEROOT 
//		  if (up_unit == 0)
//		    {
//		      proc_count = 1;
//		    }
//#endif
		}
	      
	      partial_visited_fluid_sites = 0;
	      
	      n = -1;
	      
	      for (i = 0; i < blocks_x; i++)
		{
		  for (j = 0; j < blocks_y; j++)
		    {
		      for (k = 0; k < blocks_z; k++)
			{
			  if (*(proc_id_p = &net->proc_id[ ++n ]) != marker) continue;
			  
			  if (unit_level == 1)
			    {
//#ifndef FREEROOT
			      machine_id = proc_count - net->procs;
			      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites *
							       (double)net->procs_per_machine[ machine_id ] / net->procs);
//#else
//			      machine_id = proc_count - (net->procs - 1);
//			      fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites *
//							       (double)net->procs_per_machine[ machine_id ] / (net->procs - 1));
//#endif
			    }
			  *proc_id_p = proc_count;
			  
			  if (net->id == proc_count)
			    {
			      if (vis->clusters == clusters_max)
				{
				  clusters_max <<= 1;
				  vis->cluster = (Cluster *)realloc(vis->cluster, sizeof(Cluster) * clusters_max);
				}
			      cluster_p = &vis->cluster[ vis->clusters ];
			      ++vis->clusters;
			      
			      cluster_p->block_min[0] = i;
			      cluster_p->block_min[1] = j;
			      cluster_p->block_min[2] = k;
			      
			      cluster_block_max_i = i;
			      cluster_block_max_j = j;
			      cluster_block_max_k = k;
			    }
			  partial_visited_fluid_sites += lbm->fluid_sites_per_block[ n ];
			  
			  blocks_a = 1;
			  block_location_a_p = &block_location_a[ 0 ];
			  block_location_a_p->i = i;
			  block_location_a_p->j = j;
			  block_location_a_p->k = k;
			  
			  are_fluid_sites_incrementing = 1;
			  
			  while (partial_visited_fluid_sites < fluid_sites_per_unit &&
				 are_fluid_sites_incrementing)
			    {
			      blocks_b = 0;
			      are_fluid_sites_incrementing = 0;
			      
			      for (index_a = 0;
				   index_a < blocks_a && partial_visited_fluid_sites < fluid_sites_per_unit;
				   index_a++)
				{
				  block_location_a_p = &block_location_a[ index_a ];
				  
				  for (l = 0;
				       l < 26 && partial_visited_fluid_sites < fluid_sites_per_unit;
				       l++)
				    {
				      neigh_i = block_location_a_p->i + n_x[ l ];
				      neigh_j = block_location_a_p->j + n_y[ l ];
				      neigh_k = block_location_a_p->k + n_z[ l ];
				      
				      if (neigh_i == -1 || neigh_i == blocks_x) continue;
				      if (neigh_j == -1 || neigh_j == blocks_y) continue;
				      if (neigh_k == -1 || neigh_k == blocks_z) continue;
				      
				      block_id = (neigh_i * blocks_y + neigh_j) * blocks_z + neigh_k;
				      
				      if (*(proc_id_p = &net->proc_id[ block_id ]) != marker)
					{
					  continue;
					}
				      *proc_id_p = proc_count;
				      
				      partial_visited_fluid_sites += lbm->fluid_sites_per_block[ block_id ];
				      
				      are_fluid_sites_incrementing = 1;
				      
				      if (blocks_b == blocks_buffer_size)
					{
					  blocks_buffer_size *= 2;
					  block_location_a = (BlockLocation *)realloc(block_location_a,
										      sizeof(BlockLocation) * blocks_buffer_size);
					  block_location_b = (BlockLocation *)realloc(block_location_b,
										      sizeof(BlockLocation) * blocks_buffer_size); 
					}
				      block_location_b_p = &block_location_b[ blocks_b ];
				      block_location_b_p->i = neigh_i;
				      block_location_b_p->j = neigh_j;
				      block_location_b_p->k = neigh_k;
				      ++blocks_b;
				      
				      if (proc_count == net->id)
					{
					  cluster_p->block_min[0] = min((int)cluster_p->block_min[0], neigh_i);
					  cluster_p->block_min[1] = min((int)cluster_p->block_min[1], neigh_j);
					  cluster_p->block_min[2] = min((int)cluster_p->block_min[2], neigh_k);
					  
					  cluster_block_max_i = max((int)cluster_block_max_i, neigh_i);
					  cluster_block_max_j = max((int)cluster_block_max_j, neigh_j);
					  cluster_block_max_k = max((int)cluster_block_max_k, neigh_k);
					}
				    }
				}
			      
			      if (are_fluid_sites_incrementing)
				{
				  block_location_a_p = block_location_a;
				  block_location_a = block_location_b;
				  block_location_b = block_location_a_p;
				  blocks_a = blocks_b;
				}
			    }
			  if (net->id == proc_count)
			    {
			      cluster_p->x[0] = cluster_p->block_min[0] * block_size - 0.5F * sites_x;
			      cluster_p->x[1] = cluster_p->block_min[1] * block_size - 0.5F * sites_y;
			      cluster_p->x[2] = cluster_p->block_min[2] * block_size - 0.5F * sites_z;
			      
			      cluster_p->blocks_x = 1 + cluster_block_max_i - cluster_p->block_min[0];
			      cluster_p->blocks_y = 1 + cluster_block_max_j - cluster_p->block_min[1];
			      cluster_p->blocks_z = 1 + cluster_block_max_k - cluster_p->block_min[2];
			    }
			  if (partial_visited_fluid_sites >= fluid_sites_per_unit)
			    {
			      ++proc_count;
			      partial_visited_fluid_sites = 0;
			    }
			}
		    }
		}
	    }
	}
    }
  free(block_location_b);
  free(block_location_a);
  
  net->dd_time = myClock () - seconds;
  seconds = myClock ();
  
  
  net->data_block = (DataBlock *)malloc(sizeof(DataBlock) * blocks);
  
  for (n = 0; n < blocks; n++)
    {
      net->data_block[ n ].site_data = NULL;
    }
  
  net->map_block = (DataBlock *)malloc(sizeof(DataBlock) * blocks);
  
  for (n = 0; n < blocks; n++)
    {
      net->map_block[ n ].site_data = NULL;
    }
  
  // the non-void blocks of the reference processor ("net->id") and
  // their neighbouring ones are allocated
  
  n = -1;
  
  for (i = 0; i < blocks_x; i++)
    {
      for (j = 0; j < blocks_y; j++)
	{
	  for (k = 0; k < blocks_z; k++)
	    {
	      proc_id = net->proc_id[ ++n ];
	      
	      if (proc_id != net->id) continue;
	      
	      if (net->data_block[ n ].site_data == NULL)
		{
		  net->data_block[ n ].site_data =
		    (unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
		}
	      
	      for (l = 0; l < 26; l++)
		{
		  neigh_i = i + n_x[ l ];
		  neigh_j = j + n_y[ l ];
		  neigh_k = k + n_z[ l ];
		  
		  if (neigh_i == -1 || neigh_i == blocks_x) continue;
		  if (neigh_j == -1 || neigh_j == blocks_y) continue;
		  if (neigh_k == -1 || neigh_k == blocks_z) continue;
		  
		  block_id = (neigh_i * blocks_y + neigh_j) * blocks_z + neigh_k;
		  
		  if (lbm->fluid_sites_per_block[ block_id ] == 0)
		    {
		      continue;
		    }
		  if (net->data_block[ block_id ].site_data == NULL)
		    {
		      net->data_block[ block_id ].site_data =
			(unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
		    }
		}
	    }
	}
    }
  
  // the configurational file is read and the information of non-void
  // sites is stored;
  
  site_data = (unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
  
  net->fr_time = myClock ();
  
  FILE *system_config = fopen (lbm->system_file_name, "r");
  
  xdrstdio_create (&xdr_config, system_config, XDR_DECODE);
  
  xdr_double (&xdr_config, &lbm->lattice_to_system);
  xdr_int    (&xdr_config, &blocks_x);
  xdr_int    (&xdr_config, &blocks_y);
  xdr_int    (&xdr_config, &blocks_z);
  xdr_int    (&xdr_config, &block_size);
  
  net->my_sites = 0;
  
  for (n = 0; n < blocks; n++)
    {
      xdr_int (&xdr_config, &flag);
      
      if (flag == 0) continue;
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  xdr_u_int (&xdr_config, &site_data[ m ]);
	}
      data_block_p = &net->data_block[ n ];
      
      if (data_block_p->site_data == NULL)
	{
	  continue;
	}
      proc_id = net->proc_id[ n ];
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  data_block_p->site_data[ m ] = site_data[ m ];
	  
	  if (proc_id == net->id && (site_data[ m ] & SITE_TYPE_MASK) != SOLID_TYPE)
	    {
	      ++net->my_sites;
	    }
	}
    }
  fclose (system_config);
  
  net->fr_time = myClock () - net->fr_time;
  
  free(site_data);
  
  // a map between the two-level data representation and the 1D
  // compact one is created here
  
  site_data = (unsigned int *)malloc(sizeof(unsigned int) * net->my_sites);
  
  my_sites = 0;
  
  for (n = 0; n < blocks; n++)
    {
      data_block_p = &net->data_block[ n ];
      
      if (data_block_p->site_data == NULL)
	{
	  continue;
	}
      map_block_p = &net->map_block[ n ];
      map_block_p->site_data = (unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
      
      proc_id = net->proc_id[ n ];
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  if (proc_id == net->id)
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
	    }
	  else
	    {
	      if ((data_block_p->site_data[ m ] & SITE_TYPE_MASK) == SOLID_TYPE)
		{
		  map_block_p->site_data[ m ] = (1U << 31U);
		}
              else
                {
                  map_block_p->site_data[ m ] = 0U;
                }
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
  
  my_sites = 0;
  
  n = -1;
  
  for (i = 0; i < sites_x; i += block_size)
    {
      for (j = 0; j < sites_y; j += block_size)
	{
	  for (k = 0; k < sites_z; k += block_size)
	    {
	      proc_id = net->proc_id[ ++n ];
	      
	      if (proc_id != net->id) continue;
	      
	      map_block_p = &net->map_block[ n ];
	      
	      m = -1;
	      
	      for (site_i = i; site_i < i + block_size; site_i++)
		{
		  for (site_j = j; site_j < j + block_size; site_j++)
		    {
		      for (site_k = k; site_k < k + block_size; site_k++)
			{
			  if (map_block_p->site_data[ ++m ] & (1U << 31U))
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
			      
			      site_data_p = netSiteMapPointer (neigh_i, neigh_j, neigh_k, net);
			      
			      if (site_data_p == NULL || *site_data_p & (1U << 31U))
			      	{
			      	  continue;
			      	}
			      neigh_proc_id = *netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
			      
			      if (neigh_proc_id != net->id)
				{
				  is_inner_site = 0;
				  is_inter_site = 1;
				}
			      else
				{
				  continue;
				}
			      for (mm = 0, flag = 1; mm < net->neigh_procs && flag; mm++)
				{
				  neigh_proc_p = &net->neigh_proc[ mm ];
				  
				  if (neigh_proc_id == neigh_proc_p->id)
				    {
				      flag = 0;
				      ++neigh_proc_p->fs;
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
				  neigh_proc_p->id = neigh_proc_id;
				  ++neigh_proc_p->fs;
				  ++net->neigh_procs;
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
  
  for (m = 0; m < COLLISION_TYPES; m++)
    {
      net->my_inter_collisions_sse[m] = (net->my_inter_collisions[m]/SIMD_SIZE)*SIMD_SIZE;
      net->my_inner_collisions_sse[m] = (net->my_inner_collisions[m]/SIMD_SIZE)*SIMD_SIZE;
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
  
  n = -1;
  
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
  
  // the precise interface-dependent data (interface-dependent fluid
  // site locations and identifiers of the distribution functions
  // streamed between different partitions) are collected and the
  // buffers needed for the communications are set from here
  
  net->shared_fs = 0;
  
  for (n = 0; n < net->neigh_procs; n++)
    {
      net->neigh_proc[ n ].f_data    = &f_data[ net->shared_fs<<2 ];
      net->neigh_proc[ n ].f_to_send = &f_to_send[ net->shared_fs ];
      net->neigh_proc[ n ].f_to_recv = &f_to_recv[ net->shared_fs ];
      net->neigh_proc[ n ].f_send_id = &f_send_id[ net->shared_fs ];
      net->neigh_proc[ n ].f_recv_iv = &f_recv_iv[ net->shared_fs ];
      
      m = net->neigh_proc[ n ].fs;
      net->neigh_proc[ n ].d_to_send_p = (double **)malloc(sizeof(double *) * m);
      
      net->shared_fs += net->neigh_proc[ n ].fs;
      net->neigh_proc[ n ].fs = 0;
    }
  if (net->shared_fs >= SHARED_DISTRIBUTIONS_MAX)
    {
      printf (" too many shared distributions\n");
      printf (" the execution is terminated\n");
#ifndef NOMPI
      net->err = MPI_Abort (MPI_COMM_WORLD, 1);
#else
      exit(1);
#endif
    }
  
  f_old = (double *)malloc(sizeof(double) * (net->my_sites * 15 + 1));
  f_new = (double *)malloc(sizeof(double) * (net->my_sites * 15 + 1));
  
  if (net->my_sites > 0)
    {
      f_id = (int *)malloc(sizeof(int) * (net->my_sites * 15));
      
      net->site_data = (unsigned int *)malloc(sizeof(unsigned int) * net->my_sites);
    }
  d = (double *)malloc(sizeof(double) * (net->my_sites + 1));
  
  nd_p = (double **)malloc(sizeof(double *) * (net->my_sites * 14 + 1));
  
#ifndef NOMPI
  net->err = MPI_Gather (&net->my_sites, 1, MPI_INT, proc_sites, 1, MPI_INT, 
			 0, MPI_COMM_WORLD);
#endif
  
  
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
	      proc_id = net->proc_id[ ++n ];
	      
	      if (proc_id != net->id) continue;
	      
	      map_block_p  = &net->map_block[ n ];
	      
	      m = -1;
	      
	      for (site_i = i; site_i < i + block_size; site_i++)
		{
		  for (site_j = j; site_j < j + block_size; site_j++)
		    {
		      for (site_k = k; site_k < k + block_size; site_k++)
			{
			  site_map = map_block_p->site_data[ ++m ];
			  
			  if (site_map & (1U << 31U))
			    {
			      continue;
			    }
			  f_id[ site_map*15+0 ] = site_map * 15 + 0;
			  
			  for (l = 1; l < 15; l++)
			    {
			      neigh_i = site_i + e_x[ l ];
			      neigh_j = site_j + e_y[ l ];
			      neigh_k = site_k + e_z[ l ];
			      
			      f_id[ site_map*15+l   ] = net->my_sites * 15;
			      nd_p[ site_map*14+l-1 ] = &d[ net->my_sites ];
			      
			      site_data_p = netSiteMapPointer (neigh_i, neigh_j, neigh_k, net);
			      
			      if (site_data_p == NULL || *site_data_p & (1U << 31U))
				{
				  continue;
				}
			      neigh_proc_id = *netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
			      
			      if (neigh_proc_id == net->id)
				{
				  f_id[ site_map*15+l   ] = *site_data_p * 15 + l;
				  nd_p[ site_map*14+l-1 ] = &d[ *site_data_p ];
				  continue;
				}
			      
			      neigh_proc_index = net->from_proc_id_to_neigh_proc_index[ neigh_proc_id ];
			      
			      neigh_proc_p = &net->neigh_proc[ neigh_proc_index ];
			      
			      f_data_p = &neigh_proc_p->f_data[ neigh_proc_p->fs<<2 ];
			      f_data_p[ 0 ] = site_i;
			      f_data_p[ 1 ] = site_j;
			      f_data_p[ 2 ] = site_k;
			      f_data_p[ 3 ] = l;
			      ++neigh_proc_p->fs;
			    }
			  if (!lbm->is_checkpoint)
			    {
			      d[ site_map ] = (double)lbm->block_density[ lbm->block_map[n] ];
			    }
			  net->site_data[ site_map ] = site_data[ my_sites ];
			  ++my_sites;
			}
		    }
		}
	    }
	}
    }
  free(site_data);
  
  free(lbm->block_map);
  free(lbm->block_density);
  
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
      net->req[ m ] = (MPI_Request *)malloc(sizeof(MPI_Request) * (2 * net->procs * net->procs));
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
				&net->req[ 0 ][ net->id * net->procs + m ]);
#endif
  	}
      else
  	{
#ifndef NOMPI
  	  net->err = MPI_Irecv (&neigh_proc_p->f_data[ 0 ],
  				neigh_proc_p->fs * 4, MPI_SHORT,
  				neigh_proc_p->id, 10, MPI_COMM_WORLD,
  				&net->req[ 0 ][ (net->id + net->procs) * net->procs + m ]);
#endif
  	}
    }
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      if (neigh_proc_p->id > net->id)
  	{
#ifndef NOMPI
  	  net->err = MPI_Wait (&net->req[ 0 ][ net->id * net->procs + m ], net->status);
#endif
  	}
      else
  	{
#ifndef NOMPI
  	  net->err = MPI_Wait (&net->req[ 0 ][ (net->id + net->procs) * net->procs + m ], net->status);
#endif
  	  
  	  for (n = 0; n < (neigh_proc_p->fs<<2); n += 4)
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
  for (m = 0; m < net->neigh_procs; m++)
    {
      neigh_proc_p = &net->neigh_proc[ m ];
      
      for (n = 0; n < neigh_proc_p->fs; n++)
	{
	  f_data_p = &neigh_proc_p->f_data[ n<<2 ];
	  i = f_data_p[ 0 ];
	  j = f_data_p[ 1 ];
	  k = f_data_p[ 2 ];
	  l = f_data_p[ 3 ];
	  
	  site_map = *netSiteMapPointer (i, j, k, net);
	  neigh_proc_p->f_send_id[ n ] = site_map * 15 + l;
	  neigh_proc_p->f_recv_iv[ n ] = site_map * 15 + inv_dir[ l ];
	  
	  neigh_proc_p->d_to_send_p[ n ] = &d[ site_map ];
	  nd_p[ site_map*14+l-1 ] = &neigh_proc_p->f_to_recv[ n ];
	}
    }
  
  for (n = 0; n < blocks; n++)
    {
      if (net->proc_id[ n ] != net->id)
	{
	  if (net->map_block[ n ].site_data != NULL)
	    {
	      free(net->map_block[ n ].site_data);
	      net->map_block[ n ].site_data = NULL;
	    }
	}
    }
  
  if (net->my_sites > 0)
    {
      vel = (double *)malloc(sizeof(double) * 3 * net->my_sites);
      
      for (i = 0; i < net->my_sites; i++)
	{
	  vel[ 3*i+0 ] = 1.e+30;
	  vel[ 3*i+1 ] = 1.e+30;
	  vel[ 3*i+2 ] = 1.e+30;
	}
      
      flow_field = (float *)malloc(sizeof(float) * 3 * net->my_sites);
    }
  
  net->bm_time = myClock () - seconds - net->fr_time;
}


void netEnd (Net *net, Vis *vis)
{
  // the allocated data are freed with this function
  
  int i;
  
  
  free(net->from_proc_id_to_neigh_proc_index);
  net->from_proc_id_to_neigh_proc_index = NULL;
  
  if (net->my_inner_sites + net->my_inter_sites > 0)
    {
      free(flow_field);
    }
  
  free(net->proc_id);
  net->proc_id = NULL;
  
  for (i = 0; i < blocks; i++)
    {
      if (net->map_block[ i ].site_data != NULL)
	{
	  free(net->map_block[ i ].site_data);
	  net->map_block[ i ].site_data = NULL;
	}
    }
  free(net->map_block);
  net->map_block = NULL;
  
  if (net->my_inner_sites + net->my_inter_sites > 0)
    {
      free(vel);
      vel = NULL;
      
      free(net->site_data);
      net->site_data = NULL;
      
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
  
  free(net->procs_per_machine);
  free(net->machine_id);
}

