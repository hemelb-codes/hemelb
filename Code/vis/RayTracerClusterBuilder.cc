#include "debug/Debugger.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "vis/RayTracer.h"
#include "vis/RayTracerClusterBuilder.h"

#include "log/Logger.h"

namespace hemelb
{
  namespace vis
  {
    RayTracerClusterBuilder:: RayTracerClusterBuilder
    (const geometry::LatticeData*& iLatDat,
     std::vector<Cluster> & i_clusters,
     float **& i_cluster_voxel,
     float ***& i_cluster_flow_field
     ) :
      m_clusters(i_clusters),
      m_cluster_voxel(i_cluster_voxel),
      m_cluster_flow_field(i_cluster_flow_field),
      mLatDat(iLatDat)
    {

      //Each block has a cluster id, initially set to -1
      m_cluster_id_of_block = new short int[mLatDat->GetBlockCount()];
      for (unsigned int n = 0;
	   n < mLatDat->GetBlockCount(); 
	   n++)
	{
	  m_cluster_id_of_block[n] = -1;
	}
    }
 
    RayTracerClusterBuilder::~RayTracerClusterBuilder()
    {
        delete[] m_cluster_id_of_block;
    }
    
    void RayTracerClusterBuilder::BuildClusters()
    {
      LocateClusters();
  
      m_cluster_voxel = new float *[mLatDat->GetLocalFluidSiteCount() * VIS_FIELDS];

      m_cluster_flow_field = new float **[m_clusters.size()];

      // For every cluster
      for (unsigned int lThisClusterId = 0; lThisClusterId < m_clusters.size(); lThisClusterId++)
	{
	   ProcessCluster(lThisClusterId);	  
	}


    }

    void RayTracerClusterBuilder::LocateClusters()
    {
      //Initially no blocks have been visited
      m_is_block_visited = new bool[mLatDat->GetBlockCount()];
      for (unsigned int n = 0; n < mLatDat->GetBlockCount(); n++)
	{
	  m_is_block_visited[n] = false;
	}


      // Run through all unvisited blocks finding clusters
      // Blocks are visited within the method if assigned to
      // a cluster
      site_t n = -1;
      for (site_t i = 0; i < mLatDat->GetXBlockCount(); i++)
	{
	  for (site_t j = 0; j < mLatDat->GetYBlockCount(); j++)
	    {
	      for (site_t k = 0; k < mLatDat->GetZBlockCount(); k++)
		{
		  n++;
		  if(!m_is_block_visited[n])
		    {
		      FindNewClusterStartingHere(
						 n,
						 Location(i, j, k));
		    }
		}
	    }
	}

      delete[] m_is_block_visited;
    }

    void RayTracerClusterBuilder::FindNewClusterStartingHere(site_t n, Location i_start_location)
    { 
      //Obtain the block
      geometry::LatticeData::BlockData * lBlock = mLatDat->GetBlock(n);

      //Return if there's nothing in the block assigned to the processor
      if (!SitesAssignedToLocalProcessorInBlock(lBlock))
	{
	  return;
	}

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("get0"); 

      m_current_cluster_min = i_start_location;
      m_current_cluster_max = i_start_location;

      m_sites_to_visit.push(i_start_location);

      //Loop over the cluster via neighbours until
      //all blocks have been visited
      while (!m_sites_to_visit.empty())
	{
	  //Get the location
	  m_current_location = m_sites_to_visit.top();
	  m_sites_to_visit.pop();

	  site_t block_id = mLatDat
	    ->GetBlockIdFromBlockCoords(    				     m_current_location.i,
									     m_current_location.j,				             m_current_location.k);

	  //Mark visited
	  m_is_block_visited[block_id] = true;

	  //Update the maximum and minimum x,y and z
	  //range of the cluster
	  UpdateMaxAndMin();

	  //Update the cluster id of the given block
	  m_cluster_id_of_block[block_id] = (short int) m_clusters.size();

	  //Check all the neighbouring blocks to see if they need visiting. Add them to the stack.
	  CheckNeighbouringBlocks();
	}
      

      
    }

    void RayTracerClusterBuilder::CheckNeighbouringBlocks()
    {
      //n_x n_y and n_z refer to all neigbouring sites, 
      //ie all combination of -1, 0  and 1 except (0,0,0)
      const int n_x[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, +0, +0, +0, +0, +0, +0, +0, +0, +1,
			  +1, +1, +1, +1, +1, +1, +1, +1 };
      const int n_y[] = { -1, -1, -1, +0, +0, +0, +1, +1, +1, -1, -1, -1, +0, +0, +1, +1, +1, -1,
			  -1, -1, +0, +0, +0, +1, +1, +1 };
      const int n_z[] = { -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +1, -1, +0, +1, -1,
			  +0, +1, -1, +0, +1, -1, +0, +1 };

      // Loop over all neighbouring blocks
      for (int l = 0; l < 26; l++)
	{
	  Location neigh_block = 
	    Location(
		     m_current_location.i + n_x[l],
		     m_current_location.j + n_y[l],
		     m_current_location.k + n_z[l]
		     );
		 
	  if (BlockValidAndNeedsVisiting(neigh_block))
	    {
	      m_sites_to_visit.push(neigh_block);
	    }
	}

      StoreCluster();
    }

    bool RayTracerClusterBuilder::BlockValidAndNeedsVisiting(Location block)
    {
           
      if (!mLatDat->IsValidBlockSite( 
				     block.i,
				     block.j,
				     block.k))
	{
	  return false;
	}  

      site_t block_id = mLatDat->GetBlockIdFromBlockCoords(
							   block.i,
							   block.j,
							   block.k);

      if (m_is_block_visited[block_id])
	{
	  return false;
	}
      
      geometry::LatticeData::BlockData * lBlock = mLatDat->GetBlock(block_id);		  
      return SitesAssignedToLocalProcessorInBlock(lBlock);
    }

    void RayTracerClusterBuilder::UpdateMaxAndMin()
    {
      m_current_cluster_min.i = 
	util::NumericalFunctions::min(m_current_cluster_min.i, m_current_location.i);
      m_current_cluster_min.j = 
	util::NumericalFunctions::min(m_current_cluster_min.j, m_current_location.j);
      m_current_cluster_min.k = 
	util::NumericalFunctions::min(m_current_cluster_min.k, m_current_location.k);

      m_current_cluster_max.i = 
	util::NumericalFunctions::max(m_current_cluster_max.i, m_current_location.i);
      m_current_cluster_max.j = 
	util::NumericalFunctions::max(m_current_cluster_max.j, m_current_location.j);
      m_current_cluster_max.k = 
	util::NumericalFunctions::max(m_current_cluster_max.k, m_current_location.k);
    }

    bool RayTracerClusterBuilder::
    SitesAssignedToLocalProcessorInBlock
    (geometry::LatticeData::BlockData * iBlock)
    {
      if (iBlock->ProcessorRankForEachBlockSite == NULL)
	{
	  return false;
	}

      for (unsigned int m = 0;
	   m < mLatDat->GetSitesPerBlockVolumeUnit(); 
	   m++)
	{
	  if (topology::NetworkTopology::Instance()->GetLocalRank() == 
	      iBlock->ProcessorRankForEachBlockSite[m])
	    {
	      return true;
	    }
	}
      return false;
    }

    void RayTracerClusterBuilder::StoreCluster()
    {
      Cluster lNewCluster;
      lNewCluster.x[0] = (float) (m_current_cluster_min.i * mLatDat->GetBlockSize()) - 0.5F
	* (float) mLatDat->GetXSiteCount();
      lNewCluster.x[1] = (float) (m_current_cluster_min.j * mLatDat->GetBlockSize()) - 0.5F
	* (float) mLatDat->GetYSiteCount();
      lNewCluster.x[2] = (float) (m_current_cluster_min.k * mLatDat->GetBlockSize()) - 0.5F
	* (float) mLatDat->GetZSiteCount();

      lNewCluster.blocks_x = (unsigned short) (1 + m_current_cluster_max.i - m_current_cluster_min.i);
      lNewCluster.blocks_y = (unsigned short) (1 + m_current_cluster_max.j - m_current_cluster_min.j);
      lNewCluster.blocks_z = (unsigned short) (1 + m_current_cluster_max.k - m_current_cluster_min.k);

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("found cluster %i, %i, %i, %i, %i, %i ", 
									  lNewCluster.x[0],
									  lNewCluster.x[1],
									  lNewCluster.x[2],
									  lNewCluster.blocks_x,
									  lNewCluster.blocks_y,
									  lNewCluster.blocks_z
									  );     

      m_clusters.push_back(lNewCluster);
      m_cluster_block_mins.push_back(m_current_cluster_min);
    }


    void RayTracerClusterBuilder::ProcessCluster(unsigned int i_cluster_id)
    {
      Cluster* lCluster = &m_clusters[i_cluster_id];

      m_cluster_flow_field[i_cluster_id] = new float *[lCluster->blocks_x * 
						       lCluster->blocks_y * 
						       lCluster->blocks_z]; 
      
      site_t voxel_min[3], voxel_max[3];
	  for (int l = 0; l < 3; l++)
	    {
	      voxel_min[l] = std::numeric_limits<site_t>::max();
	      voxel_max[l] = std::numeric_limits<site_t>::min();
	    }

	  int n = -1;

	  for (site_t i = 0; i < lCluster->blocks_x; i++)
	    {
	      for (site_t j = 0; j < lCluster->blocks_y; j++)
		{
		  for (site_t k = 0; k < lCluster->blocks_z; k++)
		    {
		      ++n;

		      site_t block_coord[3];
		      Location* mins = &m_cluster_block_mins[i_cluster_id];
		      block_coord[0] = (i + mins->i) * mLatDat->GetBlockSize();
		      block_coord[1] = (j + mins->j) * mLatDat->GetBlockSize();
		      block_coord[2] = (k + mins->k) * mLatDat->GetBlockSize();

		      site_t block_id = ( (i + mins->i) * mLatDat->GetYBlockCount() + (j + mins->j))
			* mLatDat->GetZBlockCount() + (k + mins->k);

		      m_cluster_flow_field[i_cluster_id][n] = NULL;

		      if (m_cluster_id_of_block[block_id] != (short int) i_cluster_id)
			{
			  continue;
			}

		      geometry::LatticeData::BlockData * lBlock = mLatDat->GetBlock(block_id);

		      m_cluster_flow_field[i_cluster_id][n]
			= new float[mLatDat->GetSitesPerBlockVolumeUnit() * VIS_FIELDS];

		      int m = -1;

		      site_t siteLocOnBlock[3];
		      for (siteLocOnBlock[0] = 0; siteLocOnBlock[0] < mLatDat->GetBlockSize(); siteLocOnBlock[0]++)
			for (siteLocOnBlock[1] = 0; siteLocOnBlock[1] < mLatDat->GetBlockSize(); siteLocOnBlock[1]++)
			  for (siteLocOnBlock[2] = 0; siteLocOnBlock[2] < mLatDat->GetBlockSize(); siteLocOnBlock[2]++)
			    {
			      unsigned int my_site_id;
			      my_site_id = lBlock->site_data[++m];

			      if (my_site_id & BIG_NUMBER3)
				{
				  for (unsigned int l = 0; l < VIS_FIELDS; l++)
				    m_cluster_flow_field[i_cluster_id][n][m * VIS_FIELDS + l] = -1.0F;

				  continue;
				}

			      for (unsigned int l = 0; l < VIS_FIELDS; l++)
				{
				  m_cluster_flow_field[i_cluster_id][n][m * VIS_FIELDS + l] = 1.0F;
				}

			      for (unsigned int l = 0; l < VIS_FIELDS; l++)
				{
				  m_cluster_voxel[my_site_id * VIS_FIELDS + l]
				    = &m_cluster_flow_field[i_cluster_id][n][m * VIS_FIELDS + l];
				}

			      for (int l = 0; l < 3; l++)
				{
				  voxel_min[l] = util::NumericalFunctions::min(voxel_min[l], siteLocOnBlock[l]
									       + block_coord[l]);
				  voxel_max[l] = util::NumericalFunctions::max(voxel_max[l], siteLocOnBlock[l]
									       + block_coord[l]);
				}

			    } // for ii[0..2]


            

		    } // for k
		} // for j
	    } // for i

	  lCluster->minmax_x[0] = (float) voxel_min[0] - 0.5F * (float) mLatDat->GetXSiteCount();
	  lCluster->minmax_y[0] = (float) voxel_min[1] - 0.5F * (float) mLatDat->GetYSiteCount();
	  lCluster->minmax_z[0] = (float) voxel_min[2] - 0.5F * (float) mLatDat->GetZSiteCount();

	  lCluster->minmax_x[1] = (float) (voxel_max[0] + 1) - 0.5F
            * (float) mLatDat->GetXSiteCount();
	  lCluster->minmax_y[1] = (float) (voxel_max[1] + 1) - 0.5F
            * (float) mLatDat->GetYSiteCount();
	  lCluster->minmax_z[1] = (float) (voxel_max[2] + 1) - 0.5F
            * (float) mLatDat->GetZSiteCount();
	  
    }
    
  }
}
