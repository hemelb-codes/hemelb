/*! \file topology.cc
    \brief In this file the functions useful to discover the topology used and
    to create and delete the domain decomposition and the various
    buffers are defined.

    Structs are defined in config.h.  Global variables (including those
    of the struct types) are declared in config.cc.  Global coordinate
    means coordinate within the entire system, not the coordinate on
    one procesor.
*/

#include "config.h"

/*!
Low level function that finds the pointer to the rank on which a
particular site resides.  proc_id is the only member of proc_block
(member of Net) for the site at global coordinate (site_i, site_j,
site_k).  If the site is in an empty block, return NULL.
*/
int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net)
{
  int i, j, k;                               // Coordinates of a cubic block
  int ii, jj, kk;                            // Coordinates of a site within the block
  ProcBlock *proc_block_p;                   // Pointer to the block
  
  if (site_i < 0 || site_i >= sites_x ||
      site_j < 0 || site_j >= sites_y ||
      site_k < 0 || site_k >= sites_z)       // Out of the bounding box.
    return NULL;   
  
  // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
  i = site_i >> shift;
  j = site_j >> shift;
  k = site_k >> shift;
  
  proc_block_p = &net->proc_block[(i * blocks_y + j) * blocks_z + k];
  
  if (proc_block_p->proc_id == NULL)   // If an empty (solid) block is addressed
    return NULL;   
  else
    {
      // Find site coordinates within the block
      ii = site_i - (i << shift);
      jj = site_j - (j << shift);
      kk = site_k - (k << shift);
      
      // Return pointer to proc_id[site] (the only member of
      // proc_block)
      return &proc_block_p->proc_id[(((ii << shift) + jj) << shift) + kk];
    }
}

/*!
Low level function that finds a pointer to site_data (the only
member of map_block (member of net)) for the site eith global
coordinate (site_i, site_j, site_k).  If the site is in an empty
block, return NULL.
*/
unsigned int *netSiteMapPointer (int site_i, int site_j, int site_k, Net *net)
{
  int i, j, k;                               // Coordinates of a block
  int ii, jj, kk;                            // Coordinates of a site within the block
  DataBlock *map_block_p;                    // Pointer to the block

  if (site_i < 0 || site_i >= sites_x ||
      site_j < 0 || site_j >= sites_y ||
      site_k < 0 || site_k >= sites_z)       // If site is out of the bounding box.
    return NULL;   

  // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
  i = site_i >> shift;
  j = site_j >> shift;
  k = site_k >> shift;

  map_block_p = &net->map_block[(i * blocks_y + j) * blocks_z + k];

  if (map_block_p->site_data == NULL)  // if an empty (solid) block is addressed
    return NULL;   
  else
    {
      // Find site coordinates within the block
      ii = site_i - (i << shift);
      jj = site_j - (j << shift);
      kk = site_k - (k << shift);
 
      // Return pointer to site_data[site]
      return &map_block_p->site_data[(((ii << shift) + jj) << shift) + kk];
    }
}

//#undef MPICHX_TOPOLOGY_DEPTHS
#ifdef MPICHX_TOPOLOGY_DEPTHS

/*!
If one has more than one machine. The topology discovery mechanism is implemented in this function
*/
int netFindTopology (Net *net, int *depths)
{
  
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

/*!
If one has more than one machine. The topology discovery mechanism is implemented in this function
*/
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


/*!
This is called from the main function.  First function to deal with processors.
The domain partitioning technique and the management of the
buffers useful for the inter-processor communications are
implemented in this function.  The domain decomposition is based
on a graph growing partitioning technique.
*/
void netInit (LBM *lbm, Net *net)
{
  double seconds;
  
  int site_i, site_j, site_k;                // Global coordinates of a site.
  int neigh_i, neigh_j, neigh_k;             // Global coordinates of a neighbour site.
  int i, j, k;                               // Global block index.
  int l;                                     // Index for neighbours of a site.
  int m;                                     // Site index on a paricular block.
  int n;                                     // Global block index.
  int mm;                                    // Index of processors surrounding this one.
  int sites_a;                               // Sites on the edge of the cluster at the start of the
                                             // current graph growing partitioning step.
  int sites_b;                               // Sites added to the edge of the cluster during the iteration.
                                             
  int index_a;                               // Site we are starting from.
  int sites_buffer_size;
  int unvisited_fluid_sites;                 // Fluid sites not yet visited.
  int partial_visited_fluid_sites;           // Fluid sites visited on a particular processor.
  int proc_count;                            // Rank we are looking at.
  int fluid_sites_per_unit;                  // Fluid sites per rank.
  int neigh_proc_index;                      
  int my_sites;                              // Sites residing on this rank.
  int are_fluid_sites_incrementing;          // Region is not bound by solid or visited sites.
  int is_inter_site, is_inner_site;          // Whether the site is on the edge of the subdomain of a rank.
  int collision_offset[2][COLLISION_TYPES];
  int flag;                                  // Whether a neighbouring process has been listed (0) or 
                                             // not (1).
  int *proc_id_p;                            // Pointer to the rank on which a particular fluid site
                                             // resides.
  
  short int *f_data_p;
  
  unsigned int *site_data;                   // Local variable to store data for sites on this rank.
  unsigned int *site_data_p; 
  unsigned int site_map;
  
  bool *is_my_block;                         // Array to store whether any sites on a block are fluid
                                             // sites residing on this rank.
  
  //Pointers fitting into SiteLocation type struct (a struct to store coordinates of a block).
  SiteLocation *site_location_a, *site_location_b;
  SiteLocation *site_location_a_p, *site_location_b_p;
  
  DataBlock *data_block_p;                   // Pointer fitting into DataBlock type struct (a single-member struct). 
  DataBlock *map_block_p;                    // Pointer fitting into DataBlock type struct (a single-member struct).
  ProcBlock *proc_block_p;                   // Pointer fitting into ProcBlock type struct (a single-member struct).
  NeighProc *neigh_proc_p;                   // Pointer fitting into NeighProc type struct.

  // Allocations.  net->fluid sites will store actual number of fluid
  // sites per proc.  Site location will store up to 10000 of some
  // sort of coordinate.
  net->fluid_sites = (int *)malloc(sizeof(int) * net->procs);
  
  sites_buffer_size = 10000;
  site_location_a = (SiteLocation *)malloc(sizeof(SiteLocation) * sites_buffer_size);
  site_location_b = (SiteLocation *)malloc(sizeof(SiteLocation) * sites_buffer_size);
  
  net->my_sites = 0;
  
  // a fast graph growing partitioning technique which spans the data
  // set only once is implemented here; the data set is explored by
  // means of the arrays "site_location_a[]" and "site_location_b[]"

  // Find the maximum number of fluid sites per process.  If steering,
  // leave one process out.
#ifndef NO_STEER 
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
#else
  fluid_sites_per_unit = (int)ceil((double)lbm->total_fluid_sites / (double)net->procs);
  proc_count = 0;
#endif


  for (n = 0; n < net->procs; n++)
    {
      net->fluid_sites[ n ] = 0;
    }
  partial_visited_fluid_sites = 0;
  unvisited_fluid_sites = lbm->total_fluid_sites;
  
  seconds = myClock ();

  if (net_machines == 1 || net_machines == net->procs) // If one machine or one machine per proc.
    {
      n = -1;
     
      // Domain Decomposition.  Pick a site. Set it to the rank we are
      // looking at. Find its neighbours and put those on the same
      // rank, then find the next-nearest neighbours, etc. until we
      // have a completely joined region, or there are enough fluid
      // sites on the rank.  In the former case, start again at
      // another site. In the latter case, move on to the next rank.
      // Do this until all sites are assigned to a rank. There is a
      // high chance of of all sites on a rank being joined.
      
      for (i = 0; i < blocks_x; i++)
	for (j = 0; j < blocks_y; j++)
	  for (k = 0; k < blocks_z; k++)
	    {
	      // Point to a block of proc_id.  If we are in a block of solids, move on.
	      proc_block_p = &net->proc_block[ ++n ];
	      
	      if (proc_block_p->proc_id == NULL)
		{
		  continue;
		}
	      m = -1;
	      
	      for (site_i = i * block_size; site_i < i * block_size + block_size; site_i++)
		for (site_j = j * block_size; site_j < j * block_size + block_size; site_j++)
		  for (site_k = k * block_size; site_k < k * block_size + block_size; site_k++)
		    {
		      // Move on if the site is solid (proc_id = 1 << 30) or has already been assigned 
		      // to a rank (0 <= proc_id < 1 << 30).  proc_id is allocated and initialised
		      // in lbmReadConfig in io.cc.
		      if (proc_block_p->proc_id[ ++m ] != -1)
			{
			  continue;
			}
		      // We have found an unvisited fluid site to start growing the subdomain from.  
		      // Assign it to the rank and update the fluid site counters.
		      proc_block_p->proc_id[ m ] = proc_count;
		      
		      if (proc_count == net->id)
			{
			  ++net->my_sites;
			}
		      ++partial_visited_fluid_sites;			      
		      ++net->fluid_sites[ proc_count ];			      
		      sites_a = 1;
		      
		      // Record the location of this initial site.
		      site_location_a_p = &site_location_a[ 0 ];
		      site_location_a_p->i = site_i;
		      site_location_a_p->j = site_j;
		      site_location_a_p->k = site_k;
		      
		      // The subdomain can grow.
		      are_fluid_sites_incrementing = 1;
		      
		      // While the region can grow (i.e. it is not bounded by solids or visited
		      // sites), and we need more sites on this particular rank.
		      while (partial_visited_fluid_sites < fluid_sites_per_unit &&
			     are_fluid_sites_incrementing)
			{
			  sites_b = 0;
			  are_fluid_sites_incrementing = 0;
			  
			  // For sites on the edge of the domain (sites_a), deal with the neighbours.
			  for (index_a = 0;
			       index_a < sites_a && partial_visited_fluid_sites < fluid_sites_per_unit;
			       index_a++)
			    {
			      site_location_a_p = &site_location_a[ index_a ];
			      
			      for (l = 1;
				   l < 15 && partial_visited_fluid_sites < fluid_sites_per_unit;
				   l++)
				{
				  // Record neighbour location.
				  neigh_i = site_location_a_p->i + e_x[ l ];
				  neigh_j = site_location_a_p->j + e_y[ l ];
				  neigh_k = site_location_a_p->k + e_z[ l ];
				  
				  // Move on if neighbour is outside the bounding box.
				  if (neigh_i == -1 || neigh_i == sites_x) continue;
				  if (neigh_j == -1 || neigh_j == sites_y) continue;
				  if (neigh_k == -1 || neigh_k == sites_z) continue;
				  
				  // Move on if the neighbour is in a block of solids (in which case
				  // the pointer to proc_id is NULL) or it is solid or has already
				  // been assigned to a rank (in which case proc_id != -1).  proc_id
				  // was initialized in lbmReadConfig in io.cc.
				  proc_id_p = netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
				  
				  if (proc_id_p == NULL || *proc_id_p != -1)
				    {
				      continue;
				    }
				  // Set the rank for a neighbour and update the fluid site counters.
				  *proc_id_p = proc_count;					  
				  ++partial_visited_fluid_sites;				  
				  ++net->fluid_sites[ proc_count ];
				  
				  if (proc_count == net->id)
				    {
				      ++net->my_sites;
				    }
				  // Neighbour was found, so the region can grow.
				  are_fluid_sites_incrementing = 1;
				  
				  // If the new layer of neighbours is too large, allocate more 
				  // memory.
				  if (sites_b == sites_buffer_size)
				    {
				      sites_buffer_size *= 2;
				      site_location_a = (SiteLocation *)realloc(site_location_a,
										sizeof(SiteLocation) * sites_buffer_size);
				      site_location_b = (SiteLocation *)realloc(site_location_b,
										sizeof(SiteLocation) * sites_buffer_size); 
				    }
				  
				  // Record the location of the neighbour.
				  site_location_b_p = &site_location_b[ sites_b ];
				  site_location_b_p->i = neigh_i;
				  site_location_b_p->j = neigh_j;
				  site_location_b_p->k = neigh_k;
				  ++sites_b;
				}
			    }
			  // When the new layer of edge sites has been found, swap the buffers for 
			  // the current and new layers of edge sites.
			  site_location_a_p = site_location_a;
			  site_location_a = site_location_b;
			  site_location_b = site_location_a_p;
			  sites_a = sites_b;
			}
		      // If we have enough sites, we have finished.
		      if (partial_visited_fluid_sites >= fluid_sites_per_unit)
			{
			  ++proc_count;
			  unvisited_fluid_sites -= partial_visited_fluid_sites;
			  fluid_sites_per_unit = (int)ceil((double)unvisited_fluid_sites / (double)(net->procs - proc_count));
			  partial_visited_fluid_sites = 0;
			}
		      // If not, we have to start growing a different region for the same rank:
		      // region expansions could get trapped.
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
		for (j = 0; j < blocks_y; j++)
		  for (k = 0; k < blocks_z; k++)
		    {
		      proc_block_p = &net->proc_block[ ++n ];
		      
		      if (proc_block_p->proc_id == NULL)
			{
			  continue;
			}
		      m = -1;
		      
		      for (site_i = i * block_size; site_i < i * block_size + block_size; site_i++)
			for (site_j = j * block_size; site_j < j * block_size + block_size; site_j++)
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
  free(site_location_b);
  free(site_location_a);
  
  net->dd_time = myClock () - seconds;
  seconds = myClock ();
  
  // A map between the two-level data representation and the 1D
  // compact one is created here.

  // Allocate an array of structures to store the fluid site identifiers for each block.
  net->map_block = (DataBlock *)malloc(sizeof(DataBlock) * blocks);
  
  for (n = 0; n < blocks; n++)
    {
      net->map_block[ n ].site_data = NULL;
    }
  
  // Local array to store site data for this rank.
  site_data = (unsigned int *)malloc(sizeof(unsigned int) * net->my_sites);
 
  // Allocate blocks.
  is_my_block = (bool *)malloc(sizeof(bool) * blocks);
  
  for (n = 0; n < blocks; n++)
    {
      is_my_block[ n ] = 0;
    }
  my_sites = 0;
  
  for (n = 0; n < blocks; n++)
    {
      // If we are in a block of solids, move to the next block.
      data_block_p = &net->data_block[ n ];
      
      if (data_block_p->site_data == NULL)
	{
	  continue;
	}
      // If we have some fluid sites, point to proc_block and map_block.
      proc_block_p = &net->proc_block[ n ];      
      map_block_p = &net->map_block[ n ];
      map_block_p->site_data = (unsigned int *)malloc(sizeof(unsigned int) * sites_in_a_block);
      
      // map_block[n].site_data is set to the fluid site identifier on this rank or (1U << 31U) if a site is solid
      // or not on this rank.  site_data is indexed by fluid site identifier and set to the site_data.
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
  
  // Free net->data_block.
  for (n = 0; n < blocks; n++)
    {
      if (net->data_block[n].site_data != NULL)
	{
	  free(net->data_block[n].site_data);
	  net->data_block[n].site_data = NULL;
	}
    }
  free(net->data_block);
  net->data_block = NULL;
  
  // If we are in a block of solids, we set map_block[n].site_data to NULL.
  for (n = 0; n < blocks; n++)
    {
      if (is_my_block[n]) continue;
      
      free(net->map_block[n].site_data);
      net->map_block[n].site_data = NULL;
      
      if (net->wall_block[n].wall_data != NULL)
	{
	  free(net->wall_block[n].wall_data);
	  net->wall_block[n].wall_data = NULL;
	}
    }
  free(is_my_block);
  
  // The numbers of inter- and intra-machine neighbouring processors,
  // interface-dependent and independent fluid sites and shared
  // distribution functions of the reference processor are calculated
  // here.  neigh_proc is a static array that is declared in config.h.
  
  net->neigh_procs = 0;
  
  for (m = 0; m < NEIGHBOUR_PROCS_MAX; m++)
    {
      net->neigh_proc[ m ].fs = 0;   // fs within NeighProc struct within the net struct.
    }
  
  net->my_inter_sites = 0;
  net->my_inner_sites = 0;
  
  for (m = 0; m < COLLISION_TYPES; m++)
    {
      net->my_inter_collisions[ m ] = 0;
      net->my_inner_collisions[ m ] = 0;
    }
  net->shared_fs = 0;   // shared fs within Net struct.
  
  my_sites = 0;
  
  n = -1;
  
  // Here, i, j and k are not the block coordinates, but the block coords * block_size.
  for (i = 0; i < sites_x; i += block_size)
    for (j = 0; j < sites_y; j += block_size)
      for (k = 0; k < sites_z; k += block_size)
	{
	  map_block_p = &net->map_block[ ++n ];
	  
	  if (map_block_p->site_data == NULL) continue;
	  
	  proc_block_p = &net->proc_block[ n ];
	  
	  m = -1;
	  
	  for (site_i = i; site_i < i + block_size; site_i++)
	    for (site_j = j; site_j < j + block_size; site_j++)
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
		      
		      // Move on if the neighbour is in a block of solids (in which case
		      // the pointer to proc_id is NULL) or it is solid (in which case proc_id ==
		      // 1 << 30) or the neighbour is also on this rank.  proc_id was initialized
		      // in lbmReadConfig in io.cc.
		      if (proc_id_p == NULL || *proc_id_p == net->id || *proc_id_p == (1 << 30))
			{
			  continue;
			}
		      is_inner_site = 0;
		      is_inter_site = 1;
		      
		      // The first time, we set mm = 0, flag = 1, but net_neigh_procs = 0, so
		      // the loop is not executed.
		      for (mm = 0, flag = 1; mm < net->neigh_procs && flag; mm++)
			{
			  // Check whether the rank for a particular neighbour has already been
			  // used for this processor.  If it has, set flag to zero.
			  neigh_proc_p = &net->neigh_proc[ mm ];
			  
			  // If proc_id is equal to a neigh_proc that has alredy been listed.
			  if (*proc_id_p == neigh_proc_p->id)
			    {
			      flag = 0;
			      ++neigh_proc_p->fs;
			      ++net->shared_fs;
			    }
			}
		      // If flag is 1, we need a new neighbouring processor.
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
			  // Store rank of neighbour in net->>neigh_proc[net->neigh_procs]
			  neigh_proc_p = &net->neigh_proc[ net->neigh_procs ];
			  neigh_proc_p->id = *proc_id_p;
			  ++neigh_proc_p->fs;
			  ++net->neigh_procs;
			  ++net->shared_fs;
			}
		    }
		  // Collision Type set here. map_block site data is renumbered according to 
		  // fluid site numbers within a particular collision type.
		  if (lbmCollisionType (site_data[ my_sites ]) == FLUID)
		    {
		      l = 0;
		    }
		  else if (lbmCollisionType (site_data[ my_sites ]) == EDGE)
		    {
		      l = 1;
		    }
		  else if (lbmCollisionType (site_data[ my_sites ]) == INLET)
		    {
		      l = 2;
		    }
		  else if (lbmCollisionType (site_data[ my_sites ]) == OUTLET)
		    {
		      l = 3;
		    }
		  else if (lbmCollisionType (site_data[ my_sites ]) == (INLET | EDGE))
		    {
		      l = 4;
		    }
		  else if (lbmCollisionType (site_data[ my_sites ]) == (OUTLET | EDGE))
		    {
		      l = 5;
		    }
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
  // Calculte the number of each type of collision.
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
      
      // If we are in a block of solids, continue.
      if (map_block_p->site_data == NULL) continue;
      
      for (m = 0; m < sites_in_a_block; m++)
	{
	  site_data_p = &map_block_p->site_data[ m ];
	  
	  // If the site is solid, continue.
	  if (*site_data_p & (1U << 31U)) continue;
	  
	  // 0th collision type for inner sites, so don't do anything.
	  if (*site_data_p < 500000000)
	    {
	      continue;
	    }
	  // Renumber the sites in map_block so that the numbers are compacted together.  We have
	  // collision offset to tell us when one collision type ends and another starts.
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
  
  // Allocate f_old and f_new according to the number of sites on the process.  The extra site
  // is there for when we would stream into a solid site during the simulation, which avoids
  // an if condition at every timestep at every boundary site.  We also allocate space for the
  // shared distribution functions.  We need twice as much space when we check the convergence
  // and the extra distribution functions are
  if (!check_conv)
    {
      f_old = (double *)malloc(sizeof(double) * (net->my_sites * 15 + 1 + net->shared_fs));
      f_new = (double *)malloc(sizeof(double) * (net->my_sites * 15 + 1 + net->shared_fs));
    }
  else
    {
      f_old = (double *)malloc(sizeof(double) * (net->my_sites * 30 + 15 + 1));
      f_new = (double *)malloc(sizeof(double) * (net->my_sites * 30 + 15 + 1));
    }
  // the precise interface-dependent data (interface-dependent fluid
  // site locations and identifiers of the distribution functions
  // streamed between different partitions) are collected and the
  // buffers needed for the communications are set from here
  
  f_data = (short int *)malloc(sizeof(short int) * 4 * net->shared_fs);
  
  if (check_conv)
    {
      f_to_send = (double *)malloc(sizeof(double) * net->shared_fs*2);
      f_to_recv = (double *)malloc(sizeof(double) * net->shared_fs*2);
      
      f_send_id = (int *)malloc(sizeof(int) * net->shared_fs);
    }

  // Allocate the index in which to put the distribution functions received from the other
  // process.
  f_recv_iv = (int *)malloc(sizeof(int) * net->shared_fs);
  
  // Reset to zero again.
  net->shared_fs = 0;
  
  for (n = 0; n < net->neigh_procs; n++)
    {
      // f_data compacted according to number of shared f_s on each process.
      // f_data will be set later.
      net->neigh_proc[ n ].f_data = &f_data[ net->shared_fs<<2 ];
      
      // Pointing to a few things, but not setting any variables.
      if (!check_conv)
	{
	  // f_head points to start of shared_fs.
	  net->neigh_proc[ n ].f_head = net->my_sites * 15 + 1 + net->shared_fs;
	}
      else
	{
	  // Points to the start of the shared_fs.  
	  net->neigh_proc[ n ].f_to_send = &f_to_send[ net->shared_fs*2 ];
	  net->neigh_proc[ n ].f_to_recv = &f_to_recv[ net->shared_fs*2 ];
	  
	  net->neigh_proc[ n ].f_send_id = &f_send_id[ net->shared_fs ];
	}
      net->neigh_proc[ n ].f_recv_iv = &f_recv_iv[ net->shared_fs ];
      
      net->shared_fs += net->neigh_proc[ n ].fs;
      net->neigh_proc[ n ].fs = 0;// This is set back to 0.
    }
  if (net->my_sites > 0)
    {
      // f_id is allocated so we know which sites to get information from.
      f_id = (int *)malloc(sizeof(int) * (net->my_sites * 15));
      
      net_site_data = (unsigned int *)malloc(sizeof(unsigned int) * net->my_sites);
      
      if (lbm_stress_type == SHEAR_STRESS)
	{
	  net_site_nor = (double *)malloc(sizeof(double) * net->my_sites*3);
	}
    }
  net->from_proc_id_to_neigh_proc_index = (short int *)malloc(sizeof(short int) * net->procs);
  
  for (m = 0; m < net->procs; m++)
    {
      net->from_proc_id_to_neigh_proc_index[ m ] = -1;
    }
  // Get neigh_proc_index from proc_id.
  for (m = 0; m < net->neigh_procs; m++)
    {
      net->from_proc_id_to_neigh_proc_index[ net->neigh_proc[m].id ] = m;
    }
  my_sites = 0;
  
  n = -1;
  
  for (i = 0; i < sites_x; i += block_size)
    for (j = 0; j < sites_y; j += block_size)
	for (k = 0; k < sites_z; k += block_size)
	  {
	    map_block_p = &net->map_block[ ++n ];
	    
	    if (map_block_p->site_data == NULL) continue;
	    
	    proc_block_p = &net->proc_block[ n ];
	    
	    m = -1;
	    
	    for (site_i = i; site_i < i + block_size; site_i++)
	      for (site_j = j; site_j < j + block_size; site_j++)
		for (site_k = k; site_k < k + block_size; site_k++)
		  {
		    // If a site is not on this process, continue.
		    if (proc_block_p->proc_id[ ++m ] != net->id)
		      {
			continue;
		      }
		    // Get site data, which is the number of the fluid site on this proc..
		    site_map = map_block_p->site_data[ m ];
		    
		    // set f_id.
		    if (!check_conv)
		      {
			f_id[ site_map*15+0 ] = site_map * 15 + 0;
		      }
		    else
		      {
			f_id[ site_map*15+0 ] = site_map * 30 + 0;
		      }
		    for (l = 1; l < 15; l++)
		      {
			// Work out positions of neighbours.
			neigh_i = site_i + e_x[ l ];
			neigh_j = site_j + e_y[ l ];
			neigh_k = site_k + e_z[ l ];
			
			// initialize f_id to the rubbish site.
			if (!check_conv)
			  {
			    f_id[ site_map*15+l ] = net->my_sites * 15;
			  }
			else
			  {
			    f_id[ site_map*15+l ] = net->my_sites * 30;
			  }
			// You know which process the neighbour is on.
			proc_id_p = netProcIdPointer (neigh_i, neigh_j, neigh_k, net);
			
			if (proc_id_p == NULL || *proc_id_p == 1 << 30)
			  {
			    continue;
			  }
			// Pointer to the neihgbour.
			site_data_p = netSiteMapPointer (neigh_i, neigh_j, neigh_k, net);
			
			// If on the same proc, set f_id of the
			// current site and direction to the
			// site and direction that it sends to.
			// If we check convergence, the data for
			// each site is split into that for the
			// current and previous cycles.
			if (*proc_id_p == net->id)
			  {
			    if (!check_conv)
			      {
				f_id[ site_map*15+l ] = *site_data_p * 15 + l;
			      }
			    else
			      {
				f_id[ site_map*15+l ] = *site_data_p * 30 + l;
			      }
			    continue;
			  }
			neigh_proc_index = net->from_proc_id_to_neigh_proc_index[ *proc_id_p ];
			
			// You have neigh proc again.
			neigh_proc_p = &net->neigh_proc[ neigh_proc_index ];
			
			// This stores some coordinates.  We
			// still need to know the site number.
			// net->neigh_proc[ n ].f_data is now
			// set as well, since this points to
			// f_data.  Every process has data for
			// its neighbours which say which sites
			// on this process are shared with the
			// neighbour.
			f_data_p = &neigh_proc_p->f_data[ neigh_proc_p->fs<<2 ];
			f_data_p[ 0 ] = site_i;
			f_data_p[ 1 ] = site_j;
			f_data_p[ 2 ] = site_k;
			f_data_p[ 3 ] = l;
			++neigh_proc_p->fs; // We recount this again.
		      }
		    // This is used in Calculate BC in IO.
		    net_site_data[ site_map ] = site_data[ my_sites ];
		    
		    if (lbm_stress_type == SHEAR_STRESS)
		      {
			net_site_nor[ site_map*3 ] = 1.0e+30;
			
			if (lbmCollisionType (net_site_data[ site_map ]) & EDGE)
			  {
			    for (l = 0; l < 3; l++)
			      net_site_nor[ site_map*3+l ] = net->wall_block[n].wall_data[m].wall_nor[l];
			  }
			else
			  {
			    net_site_nor[ site_map*3 ] = 1.0e+30;
			  }
		      }
		    ++my_sites;
		  }
	  }
  free(site_data);
  
  // point-to-point communications are performed to match data to be
  // sent to/receive from different partitions; in this way, the
  // communication of the locations of the interface-dependent fluid
  // sites and the identifiers of the distribution functions which
  // propagate to different partitions is avoided (only their values
  // will be communicated). It's here!
  
  // Allocate the request variable.
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
   
      // One way send receive.  The lower numbered procs send and the higher numbered ones receive.
      // It seems that, for each pair of processors, the lower numbered one ends up with its own
      // edge sites and directions stored and the higher numbered one ends up with those on the
      // other processor.
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
	  // Now we sort the situation so that each process has its own sites.
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
	  // Get coordinates and direction of the distribution function to be sent to another process.
	  f_data_p = &neigh_proc_p->f_data[ n*4 ];
	  i = f_data_p[ 0 ];
	  j = f_data_p[ 1 ];
	  k = f_data_p[ 2 ];
	  l = f_data_p[ 3 ];
	  
	  // Get the fluid site number of site that will send data to another process.
	  site_map = *netSiteMapPointer (i, j, k, net);
	  
	  if (!check_conv)
	    {
	      // Set f_id to the element in the send buffer that we put the updated
	      // distribution functions in.
	      f_id[ site_map * 15 + l ] = ++f_count;

	      // Set the place where we put the received distribution functions, which is
	      // f_new[number of fluid site that sends, inverse direction].
	      neigh_proc_p->f_recv_iv[ n ] = site_map * 15 + inv_dir[ l ];
	    }
	  else
	    {
	      // Set f_send_id to the element of f_old that we pull the post-collisional
	      // distributions from.  f_id will send the updated distribution functions to the
	      // rubbish site instead of automatically putting them in the send buffer.
	      neigh_proc_p->f_send_id[ n ] = site_map * 30 + l;

	      // Set the place where we put the received distribution functions, which is
	      // f_new[number of fluid site that sends, inverse direction].
	      neigh_proc_p->f_recv_iv[ n ] = site_map * 30 + inv_dir[ l ];
	    }
	}
    }
  // neigh_prc->f_data was only set as a pointer to f_data, not allocated.  In this line, we 
  // are freeing both of those.
  free(f_data);
  
  net->bm_time = myClock () - seconds;
}

/*!
Free the allocated data.
*/
void netEnd (Net *net)
{
  int i;
  
  
  free(net->from_proc_id_to_neigh_proc_index);
  net->from_proc_id_to_neigh_proc_index = NULL;
  
  free(f_recv_iv);
  
  if (check_conv)
    {
      free(f_send_id);
      free(f_to_recv);
      free(f_to_send);
    }
  
  if (lbm_stress_type == SHEAR_STRESS && net->my_sites > 0)
    {
      free(net_site_nor);
      
      for (i = 0; i < blocks; i++)
	{
	  if (net->wall_block[ i ].wall_data != NULL)
	    {
	      free(net->wall_block[ i ].wall_data);
	    }
	}
      free(net->wall_block);
    }
  for (i = 0; i < blocks; i++)
    {
      if (net->map_block[ i ].site_data != NULL)
	{
	  free(net->map_block[ i ].site_data);
	}
    }
  free(net->map_block);
  
  for (i = 0; i < blocks; i++)
    {
      if (net->proc_block[ i ].proc_id != NULL)
	{
	  free(net->proc_block[ i ].proc_id);
	}
    }
  free(net->proc_block);
  
  if (net->my_sites > 0)
    {
      free(net_site_data);
      free(f_id);
    }
  free(f_new);
  free(f_old);
  
#ifndef NOMPI
  for (i = 0; i < COMMS_LEVELS; i++)
    free(net->req[ i ]);
  free(net->req);
#endif
  
  free(net->fluid_sites);
  
  free(net->procs_per_machine);
  free(net->machine_id);
}

