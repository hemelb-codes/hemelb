#ifndef HEMELB_VIS_RAYTRACERCLUSTERBUILDER_H
#define HEMELB_VIS_RAYTRACERCLUSTERBUILDER_H

#include <vector>
#include <stack>

#include "constants.h"
#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"
#include "vis/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    class RayTracerClusterBuilder
    {
    public:
      RayTracerClusterBuilder
	(const geometry::LatticeData*& iLatDat,
	 std::vector<Cluster> & i_clusters,
	 float **& i_cluster_voxel,
	 float ***& i_cluster_flow_field
	 );
      ~RayTracerClusterBuilder();
      void BuildClusters();


    private:
      void LocateClusters();
      
      // If the site hasn't been visited, finds a new rectangular
      // cluster containing this site
      void FindNewClusterStartingHere
	(site_t n, 
	 Location i_start_location);
      void CheckNeighbouringBlocks();

      bool BlockValidAndNeedsVisiting(Location block);
      bool SitesAssignedToLocalProcessorInBlock
	(geometry::LatticeData::BlockData * iBlock);

      void UpdateClusterMaxAndMin();

      void StoreCluster();

      void ProcessCluster(unsigned int i_cluster_id);

      void UpdateFlowField
	(geometry::LatticeData::BlockData * i_block,
	 unsigned int n, unsigned int i_cluster_id, int l_site_id);

      Location GetSiteCoordinatesOfBlock
	(unsigned int i_cluster_id, Location offset);

      static void UpdateMinLocation(Location io_store_location, Location i_compare_location);

      static void UpdateMaxLocation(Location io_store_location, Location i_compare_location);

      

      bool* m_is_block_visited;

   
      //The algorithm stores a list of locations on route 
      //effectively performing a depth first search
      //on the cluster
      //The algorithm operates at block level, so 
      //all locations are in terms of blocks 

      std::stack<Location> m_sites_to_visit;
      Location m_current_location;

      Location m_current_cluster_min;
      Location m_current_cluster_max;

      //The clusters are the exception to the rule -
      //the centre is stored in terms of number of
      //sites as a float float, though the extent is
      //number of blocks.
      std::vector<Cluster> & m_clusters;
      float **& m_cluster_voxel;
      float ***& m_cluster_flow_field;

      //These minimums are stored in terms of block
      //location - essentially duplicating the 
      //cluster data
      std::vector<Location> m_cluster_block_mins; 

      //These are stored in site_coordinates
      Location m_current_min_voxel;
      Location m_current_max_voxel;

      const geometry::LatticeData*& mLatDat;
      short int *m_cluster_id_of_block;
      
      
    };
  }
}

#endif // HEMELB_VIS_RAYTRACERBUILDER_H
