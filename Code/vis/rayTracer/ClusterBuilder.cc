#include "debug/Debugger.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "vis/rayTracer/Location.h"
#include "vis/rayTracer/RayTracer.h"
#include "log/Logger.h"

namespace hemelb
{
    namespace vis
    {
	namespace raytracer
	{
	    RayTracer::ClusterBuilder::ClusterBuilder
	    (const geometry::LatticeData*& iLatDat,
	     std::vector<Cluster> & i_clusters,
	     float **& i_cluster_voxel,
	     float ***& i_cluster_flow_field
		) :
		mBlockIterator(iLatDat),
		m_clusters(i_clusters),
		m_cluster_voxel(i_cluster_voxel),
		m_cluster_flow_field(i_cluster_flow_field),
		mLatDat(iLatDat)
	    {

		//Each block has a cluster id, initially set to -1
		m_cluster_id_of_block = new short int[mLatDat->GetBlockCount()];
		for (site_t n = 0;
		     n < mLatDat->GetBlockCount(); 
		     n++)
		{
		    m_cluster_id_of_block[n] = -1;
		}
	    }
 
	    RayTracer::ClusterBuilder::~ClusterBuilder()
	    {
		delete[] m_cluster_id_of_block;
	    }
    
	    void RayTracer::ClusterBuilder::BuildClusters()
	    {
		LocateClusters();
  
		m_cluster_voxel = new float * [mLatDat->GetLocalFluidSiteCount() * VIS_FIELDS];

		m_cluster_flow_field = new float ** [m_clusters.size()];

		// For every cluster
		for (site_t lThisClusterId = 0; lThisClusterId < m_clusters.size(); lThisClusterId++)
		{
		    ProcessCluster(lThisClusterId);	  
		}

	    }

	    void RayTracer::ClusterBuilder::LocateClusters()
	    {
		// Run through all unvisited blocks finding clusters
		while (mBlockIterator.GoToNextUnvisitedBlock())
		{
		    //Continue if there's nothing in the block assigned to the processor
		    if (mBlockIterator.SitesAssignedToLocalProcessorInBlock())
		    {
			FindNewCluster();
		    }
		}
	    }	 

	    void RayTracer::ClusterBuilder::FindNewCluster()
	    { 
	      //Set the cluster minimum and maximum location to the current location
	      //These will eventually contain the bounds of the rectangular cluster
	      Location lClusterMin = mBlockIterator.GetLocation();
	      Location lClusterMax = mBlockIterator.GetLocation();
		
	      //To discover the cluster, we must continually visit the neighbours of sequential blocks 
	      //We keep a stack of all the sites that must be visited
	      std::stack<Location> lBlocksToVisit;
	      lBlocksToVisit.push(mBlockIterator.GetLocation());

	      //Loop over the cluster via neighbours until
	      //all blocks have been visited
	      while (!lBlocksToVisit.empty())
	      {
		//Get the location
		Location lCurrentLocation = lBlocksToVisit.top();
		lBlocksToVisit.pop();

		mBlockIterator.MarkBlockVisited(lCurrentLocation);
		    
		if(mBlockIterator.SitesAssignedToLocalProcessorInBlock())
		{
		  //Update the x,y and z range of the cluster
		  Location::UpdateMinLocation(lClusterMin, lCurrentLocation);
		  Location::UpdateMaxLocation(lClusterMax, lCurrentLocation);
		      
		  //TODO: pdate the maximum x, y and z site location
		  		    
		  //Update the cluster id of the given block
		  site_t block_id = mBlockIterator.CurrentNumber();
		  m_cluster_id_of_block[block_id] = (short int) m_clusters.size();
			
		  //Check all the neighbouring blocks to see if they need visiting. Add them to the stack.
		  AddNeighbouringBlocks(lCurrentLocation, lBlocksToVisit);
		}
	      }
		
	      StoreCluster(lClusterMin, lClusterMax);
	    }

	  void RayTracer::ClusterBuilder::AddNeighbouringBlocks
	  (Location iCurrentLocation, std::stack<Location> iBlocksToVisit)
	    {
		//All neibouring blocks
		const Location lNeighbour[26] =
		{
		    Location(-1, -1, -1),
		    Location(-1, -1,  0),
		    Location(-1, -1,  1),
		    
		    Location(-1,  0, -1),
		    Location(-1,  0,  0),
		    Location(-1,  0,  1),

		    Location(-1,  1, -1),
		    Location(-1,  1,  0),
		    Location(-1,  1,  1),

		    
		    Location( 0, -1, -1),
		    Location( 0, -1,  0),
		    Location( 0, -1,  1),
		    
		    Location( 0,  0, -1),
		    // 0 0 0 is same site
		    Location( 0,  0,  1),

		    Location( 0,  1, -1),
		    Location( 0,  1,  0),
		    Location( 0,  1,  1),

		    Location( 1, -1, -1),
		    Location( 1, -1,  0),
		    Location( 1, -1,  1),
		    
		    Location( 1,  0, -1),
		    Location( 1,  0,  0),
		    Location( 1,  0,  1),

		    Location( 1,  1, -1),
		    Location( 1,  1,  0),
		    Location( 1,  1,  1)
		};

		// Loop over all neighbouring blocks
		for (int l = 0; l < 26; l++)
		{
		    Location lNeigbouringBlock = iCurrentLocation + lNeighbour[l];
		 
		    if (mBlockIterator.BlockValid(lNeigbouringBlock))
		    {
			if(!mBlockIterator.BlockVisited(lNeigbouringBlock))
			{
			    iBlocksToVisit.push(lNeigbouringBlock);
			}
		    }
		}
	    }

	    void RayTracer::ClusterBuilder::StoreCluster(Location iClusterMin, Location iClusterMax)
	    {
		Cluster lNewCluster;
		lNewCluster.blockCoordinates[0] = (float) (iClusterMin.i * mLatDat->GetBlockSize()) - 0.5F
		    * (float) mLatDat->GetXSiteCount();
		lNewCluster.blockCoordinates[1] = (float) (iClusterMin.j * mLatDat->GetBlockSize()) - 0.5F
		    * (float) mLatDat->GetYSiteCount();
		lNewCluster.blockCoordinates[2] = (float) (iClusterMin.k * mLatDat->GetBlockSize()) - 0.5F
		    * (float) mLatDat->GetZSiteCount();

		lNewCluster.blocks_x = (unsigned short) (1 + iClusterMax.i - iClusterMin.i);
		lNewCluster.blocks_y = (unsigned short) (1 + iClusterMax.j - iClusterMin.j);
		lNewCluster.blocks_z = (unsigned short) (1 + iClusterMax.k - iClusterMin.k);

		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("found cluster %i, %i, %i, %i, %i, %i ", 
										    lNewCluster.blockCoordinates[0],
										    lNewCluster.blockCoordinates[1],
										    lNewCluster.blockCoordinates[2],
										    lNewCluster.blocks_x,
										    lNewCluster.blocks_y,
										    lNewCluster.blocks_z
		    );     

		m_clusters.push_back(lNewCluster);
		m_cluster_block_mins.push_back(iClusterMin);
	    }


	    void RayTracer::ClusterBuilder::ProcessCluster(site_t i_cluster_id)
	    {
		Cluster* lCluster = &m_clusters[i_cluster_id];

		m_cluster_flow_field[i_cluster_id] = new float *[lCluster->blocks_x * 
								 lCluster->blocks_y * 
								 lCluster->blocks_z]; 
      
		m_current_min_voxel = Location(std::numeric_limits<site_t>::max());
		m_current_max_voxel =  Location(std::numeric_limits<site_t>::min());

		int n = -1;

		for (site_t i = 0; i < lCluster->blocks_x; i++)
		{
		    for (site_t j = 0; j < lCluster->blocks_y; j++)
		    {
			for (site_t k = 0; k < lCluster->blocks_z; k++)
			{
			    ++n;

			    Location block_coordinates = Location(i, j, k) + m_cluster_block_mins[i_cluster_id];
			    site_t block_id = mLatDat->GetBlockIdFromBlockCoords
				(block_coordinates.i,
				 block_coordinates.j,
				 block_coordinates.k);

			    m_cluster_flow_field[i_cluster_id][n] = NULL;

			    if (m_cluster_id_of_block[block_id] != (short int) i_cluster_id)
			    {
				continue;
			    }

			    geometry::LatticeData::BlockData * lBlock = mLatDat->GetBlock(block_id);

			    m_cluster_flow_field[i_cluster_id][n]
				= new float[mLatDat->GetSitesPerBlockVolumeUnit() * VIS_FIELDS];

			    UpdateBlockMaxMinsAndFlowField(lBlock, n, i_cluster_id,block_coordinates);
			} // for k
		    } // for j
		} // for i

		lCluster->minmax_x[0] = (float) m_current_min_voxel.i - 0.5F * (float) mLatDat->GetXSiteCount();
		lCluster->minmax_y[0] = (float) m_current_min_voxel.j - 0.5F * (float) mLatDat->GetYSiteCount();
		lCluster->minmax_z[0] = (float) m_current_min_voxel.k - 0.5F * (float) mLatDat->GetZSiteCount();

		lCluster->minmax_x[1] = (float) (m_current_max_voxel.i + 1) - 0.5F
		    * (float) mLatDat->GetXSiteCount();
		lCluster->minmax_y[1] = (float) (m_current_max_voxel.j + 1) - 0.5F
		    * (float) mLatDat->GetYSiteCount();
		lCluster->minmax_z[1] = (float) (m_current_max_voxel.k + 1) - 0.5F
		    * (float) mLatDat->GetZSiteCount();
	  
	    }
    
	    void RayTracer::ClusterBuilder::UpdateBlockMaxMinsAndFlowField
	    (geometry::LatticeData::BlockData * lBlock, site_t n, site_t i_cluster_id, Location i_block_coordinates)
	    {
		int l_site_id  = -1;
    
		Location site_coordinates_of_block = i_block_coordinates * mLatDat->GetBlockSize();
		Location siteLocOnBlock;

		for (siteLocOnBlock.i = 0; siteLocOnBlock.i < mLatDat->GetBlockSize(); siteLocOnBlock.i++)
		{
		    for (siteLocOnBlock.j = 0; siteLocOnBlock.j < mLatDat->GetBlockSize(); siteLocOnBlock.j++)
		    {
			for (siteLocOnBlock.k = 0; siteLocOnBlock.k < mLatDat->GetBlockSize(); siteLocOnBlock.k++)
			{
			    ++l_site_id;
			  
			    UpdateSiteFlowField(lBlock, n, i_cluster_id, l_site_id);
			  	  
			    //Come back to l
			    UpdateMinLocation(m_current_min_voxel, siteLocOnBlock + i_block_coordinates);
			    UpdateMaxLocation(m_current_max_voxel, siteLocOnBlock + i_block_coordinates);
			}
		    } 
		}
	    }
    
	    void RayTracer::ClusterBuilder::UpdateSiteFlowField
	    ( geometry::LatticeData::BlockData * i_block,
	      site_t n, site_t i_cluster_id, int l_site_id)
	    {
		//TODO: Clean this
      
		site_t site_datum = i_block->site_data[l_site_id];
      
		//If site is solid
		if (site_datum & BIG_NUMBER3)
		{
		    for (site_t l = 0; l < VIS_FIELDS; l++)
		    {		
			m_cluster_flow_field[i_cluster_id][n][l_site_id * VIS_FIELDS + l] = -1.0F;
		    }
		}
		else {
		    for (site_t l = 0; l < VIS_FIELDS; l++)
		    {
			m_cluster_flow_field[i_cluster_id][n][l_site_id * VIS_FIELDS + l] = 1.0F;
		    }

		    for (site_t l = 0; l < VIS_FIELDS; l++)
		    {
			m_cluster_voxel[site_datum * VIS_FIELDS + l]
			    = &m_cluster_flow_field[i_cluster_id][n][l_site_id * VIS_FIELDS + l];
		    }
		}
		
	    }
    

	   

	}
    }
}
