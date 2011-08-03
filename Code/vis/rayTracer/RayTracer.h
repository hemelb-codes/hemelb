#ifndef HEMELB_VIS_RAYTRACER_H
#define HEMELB_VIS_RAYTRACER_H

#include <stack>

#include "constants.h"
#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"
#include "vis/rayTracer/Location.h"



namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      struct Cluster
      {
	float minmax_x[2], minmax_y[2], minmax_z[2];
      
	//Stores the lowest x, y and z point of the Cluster 
	float blockCoordinates[3];
      
	//Stores the size of the cluster
	unsigned short int blocks_x;
	unsigned short int blocks_y;
	unsigned short int blocks_z;
      };

      class RayTracer
      {
      public:
	// Constructor and destructor do all the usual stuff.
	RayTracer(const geometry::LatticeData* iLatDat,
		  const DomainStats* iDomainStats,
		  Screen* iScreen,
		  Viewpoint* iViewpoint,
		  VisSettings* iVisSettings);
	~RayTracer();

	// Method to update the voxel corresponding to site i with its
	// newly calculated density, velocity and stress.
	void UpdateClusterVoxel(site_t i,
				distribn_t density,
				distribn_t velocity,
				distribn_t stress);

	// Render the current state into an image.
	void Render();

      private:

	class ClusterBuilder
	{
	public:
	  ClusterBuilder
	    (const geometry::LatticeData*& iLatDat,
	     std::vector<Cluster> & i_clusters,
	     float **& i_cluster_voxel,
	     float ***& i_cluster_flow_field
	      );
	  ~ClusterBuilder();
	  void BuildClusters();
      

	private:
	  class RectangularIterator
	  {
	  public:
	    RectangularIterator();

	    Location GetLocation();
			
	    site_t CurrentNumber();

	    site_t GetNumberFromLocation(Location iLocation);

	    bool Iterate();

	    virtual site_t GetXCount() = 0;

	    virtual site_t GetYCount() = 0;

	    virtual site_t GetZCount() = 0;
			
	  protected:
	    Location mCurrentLocation;
	    site_t mCurrentNumber;
	  };

	  class SiteIterator : public RectangularIterator
	  {
	  public:
	    SiteIterator(const geometry::LatticeData * iLatticeDat, 
			 const site_t iBlockId);
			
	    virtual site_t GetXCount();

	    virtual site_t GetYCount();

	    virtual site_t GetZCount();

	  private:
	    const geometry::LatticeData * mLatticeData;
	    const site_t iBlockId;

	  };
		    
	  class BlockIterator : public RectangularIterator
	  {
	  public:
	    BlockIterator(const geometry::LatticeData * iLatDat);

	    site_t CurrentBlockNumber();
		
	    bool GoToNextUnvisitedBlock();
		
	    geometry::LatticeData::BlockData * GetBlockData();

	    SiteIterator GetSiteIterator();

	    virtual site_t GetXCount();

	    virtual site_t GetYCount();

	    virtual site_t GetZCount();

	    bool BlockValid(Location block);

	    bool BlockVisited(site_t n);

	    bool BlockVisited(Location n);

	    bool CurrentBlockVisited();

	    void MarkBlockVisited();
	    void MarkBlockVisited(site_t n);
	    void MarkBlockVisited(Location location);

	    bool SitesAssignedToLocalProcessorInBlock();

	  private:
	    bool GoToNextBlock();

	    const geometry::LatticeData * mLatticeData;

	    bool* mBlockVisited;
	  
	  };

	  void LocateClusters();
      
	  // If the site hasn't been visited, finds a new rectangular
	  // cluster containing this site
	  void FindNewCluster();
	  void AddNeighbouringBlocks(Location iCurrentLocation,
				     std::stack<Location> iBlocksToVisit);

	  bool BlockValidAndNeedsVisiting(Location block);

	  bool SitesAssignedToLocalProcessorInBlock
	    (geometry::LatticeData::BlockData * iBlock);
      
	  void StoreCluster(Location iClusterMin, Location iClusterMax);

	  void ProcessCluster(site_t i_cluster_id);

	  void UpdateBlockMaxMinsAndFlowField(geometry::LatticeData::BlockData * lBlock, site_t n,
					      site_t i_cluster_id, Location i_block_coordinates);

	  void UpdateSiteFlowField
	    (geometry::LatticeData::BlockData * i_block,
	     site_t n, site_t i_cluster_id, int l_site_id);

	  Location GetSiteCoordinatesOfBlock
	    (site_t i_cluster_id, Location offset);

	  static void UpdateMinLocation(Location io_store_location, Location i_compare_location);

	  static void UpdateMaxLocation(Location io_store_location, Location i_compare_location);

	  BlockIterator mBlockIterator;

   
	  //The algorithm stores a list of locations on route 
	  //effectively performing a depth first search
	  //on the cluster
	  //The algorithm operates at block level, so 
	  //all locations are in terms of blocks 

		    
		    
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

	struct Ray
	{
	  float Direction[3];
	  float InverseDirection[3];
	  float Length;

	  float VelocityColour[3];
	  float StressColour[3];
	  float Stress;
	  float Density;
	  float MinT;
	};

	struct AABB
	{
	  float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
	};


	void UpdateRayData(const float flow_field[3],
			   float ray_t,
			   float ray_segment,
			   Ray* bCurrentRay);

	void TraverseVoxels(const float block_min[3],
			    const float block_x[3],
			    const float voxel_flow_field[],
			    float t,
			    Ray* bCurrentRay,
			    const bool xyz_is_1[3]);

	void TraverseBlocks(const Cluster* cluster,
			    const bool xyz_Is_1[3],
			    const float ray_dx[3],
			    float **block_flow_field,
			    Ray *bCurrentRay);

	void AABBvsRay(const AABB* aabb,
		       const float inverseDirection[3],
		       const bool xyzComponentIsPositive[3],
		       float* t_near,
		       float* t_far);

	void UpdateColour(float dt, const float palette[3], float col[3]);

	void BuildClusters();

	const geometry::LatticeData* mLatDat;

	const DomainStats* mDomainStats;
	Screen* mScreen;
	Viewpoint* mViewpoint;
	VisSettings* mVisSettings;

	std::vector<Cluster> mClusters;
	float **cluster_voxel;
	float ***cluster_flow_field;

	const float mBlockSizeFloat;
	const float mBlockSizeInverse;
	const site_t block_size2, block_size3, block_size_1;
	const site_t blocks_yz;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_H
