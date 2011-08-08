#ifndef HEMELB_VIS_RAYTRACER_H
#define HEMELB_VIS_RAYTRACER_H

#include <map>
#include <stack>
#include <vector>

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

	//Stores the data about an individual voxel 
	class SiteData_t
	{
	public:
	  float Density;
	  float Velocity;
	  float Stress;

	SiteData_t(float iValue) :
	  Density(iValue), Velocity(iValue), Stress(iValue) {}
	};

	//The cluster structure stores data relating to the clusters
	//used by the RayTracaer, in an optimal format
	//Clusters are produced by the ClusterFactory
	//Caution: the data within the flow field is altered by means
	//of pointers obtained from the GetClusterVoxelDataPointer
	//method
	class Cluster
	{
	public:
	  Cluster();

	  Location<float> minSite;
	  Location<float> maxSite;

	  //Stores the lowest x, y and z block location of the Cluster 
	  //in terms of site units relative to the centre location
	  Location<float> minBlock;
        
	  //Stores the size of the cluster in terms of the number of blocks
	  unsigned short int blocksX;
	  unsigned short int blocksY;
	  unsigned short int blocksZ;

	  std::vector<std::vector<SiteData_t> > SiteData;

	};

      private:

	class ClusterBuilder
	{
	public:
	  ClusterBuilder
	    (const geometry::LatticeData*& iLatDat);
	  ~ClusterBuilder();
	  
	  void BuildClusters();
	  
	  std::vector<Cluster>& GetClusters();

	  SiteData_t* GetClusterVoxelDataPointer(site_t iSiteId);

	private:
	  //Volume tracker is used to sequentially traverse a 
	  //3D structure maintaining the index and Location 
	  //within volume
	  class VolumeTraverser
	  {
	  public:
	    VolumeTraverser();

	    Location<site_t> GetCurrentLocation();
			
	    site_t GetCurrentIndex();

	    site_t GetIndexFromLocation(Location<site_t> iLocation);

	    //Increments the index by one and update the location accordingly
	    //Returns true if successful or false if the whole volume has been
	    //traversed
	    bool TraverseOne();

	    //Virtual methods which must be defined for correct traversal
	    virtual site_t GetXCount() = 0;
	    virtual site_t GetYCount() = 0;
	    virtual site_t GetZCount() = 0;
			
	  protected:
	    Location<site_t> mCurrentLocation;
	    site_t mCurrentNumber;
	  };
	  
	  //SiteTraverse is used to traverse the sites in a speficied block
	  //within the lattice data
	  class SiteTraverser : public VolumeTraverser
	  {
	  public:
	    SiteTraverser(const geometry::LatticeData * iLatticeDat, 
			  const site_t iBlockId);
			
	    virtual site_t GetXCount();

	    virtual site_t GetYCount();

	    virtual site_t GetZCount();

	  private:
	    //Returns the block size in number of sites
	    site_t GetBlockSize();

	    const geometry::LatticeData * mLatticeData;
	    const site_t mBlockId;

	  };
	
  	  //BlockTraverser is used to traverse the blocks in a lattice sequentially.
	  //The class also contains a record of which blocks have been visited, which
	  //is neccessary for the algoritm which uses this. No locations are automatically
	  //marked visited, and methods have been created to assist with random access
	  //of the lattice data as required by the algorithm
	  class BlockTraverser : public VolumeTraverser
	  {
	  public:
	    BlockTraverser(const geometry::LatticeData * iLatDat);
	    ~BlockTraverser();

	    site_t CurrentBlockNumber();
	    
	    Location<site_t> GetSiteCoordinatesOfLowestSiteInCurrentBlock();
	
	    //Tranverses the block until the next unvisited block is reached.
	    //Returns false if the end of the Volume is reached
	    bool GoToNextUnvisitedBlock();
		
	    geometry::LatticeData::BlockData * GetCurrentBlockData();

	    geometry::LatticeData::BlockData * GetBlockDataForLocation(const Location<site_t>& iLocation);

	    site_t GetBlockSize();
	    
	    SiteTraverser GetSiteTraverserForCurrentBlock();
	    
	    SiteTraverser GetSiteTraverserForLocation(const Location<site_t>& iLocation);

	    virtual site_t GetXCount();

	    virtual site_t GetYCount();

	    virtual site_t GetZCount();

	    bool IsValidLocation(Location<site_t> block);
	    
	    bool IsCurrentBlockVisited();

	    bool IsBlockVisited(site_t n);
	    bool IsBlockVisited(Location<site_t> n);

	    void MarkCurrentBlockVisited();

	    void MarkBlockVisited(site_t n);
	    void MarkBlockVisited(Location<site_t> location);

	  private:
	    bool GoToNextBlock();

	    const geometry::LatticeData * mLatticeData;

	    bool* mBlockVisited;
	  };

	  //Locates all the clusters in the lattice structure and the 
	  void LocateClusters();
      
	  // If the site hasn't been visited, finds a new rectangular
	  // cluster containing this site
	  void FindNewCluster();

	  //Adds neighbouring blocks of the input location to the input stack
	  void AddNeighbouringBlocks(Location<site_t> iCurrentLocation,
				     std::stack<Location<site_t> >& ioBlocksToVisit);

	  //Returns true if there are sites in the given block associated with the 
	  //local processor rank 
	  bool AreSitesAssignedToLocalProcessorRankInBlock
	    (geometry::LatticeData::BlockData * iBlock);
      
	  //Adds a new cluster by taking in the required data in interget format
	  //and converting it to that used by the raytracer
	  //NB: Futher processing is required on the cluster before it can be used 
	  //by the ray tracer, which is handled by the ProcessCluster method
	  void AddCluster(Location<site_t> iClusterBlockMin,
			  Location<site_t> iClusterBlockMax,
			  Location<site_t> iClusterVoxelMin,
			  Location<site_t> iClusterVoxelMax);

	  //Adds "flow-field" data to the cluster
	  void ProcessCluster(unsigned int iClusterId);

	  void UpdateSiteData
	    (geometry::LatticeData::BlockData * lBlock, site_t n, 
	     unsigned int iClusterId, Location<site_t> i_block_coordinates);

	  void UpdateSiteDataAtSite
	    (geometry::LatticeData::BlockData * iBlock,
	     site_t n, unsigned int iClusterId, unsigned int l_site_id);

	  Location<site_t> GetSiteCoordinatesOfBlock
	    (site_t iClusterId, Location<site_t> offset);

	  SiteData_t* GetDataPointerClusterVoxelSiteId(site_t iSiteId);

	  void SetDataPointerForClusterVoxelSiteId
	    (site_t iClusterVortexSiteId, 
	     SiteData_t* iDataPointer);
	  
	  BlockTraverser mBlockTraverser;
	  
	  //Caution: the data within mClusters is altered by means
	  //of pointers obtained from the GetClusterVoxelDataPointer
	  //method. No insertion of copying must therefore take place 
	  //on mClusters once building is complete
	  std::vector<Cluster> mClusters;
	 
	  std::vector<Location<site_t> > mClusterBlockMins; 
	 
	  //This allows a cluster voxel site ID (as part of the 1D structure for)
	  //storing sites to be mapped to the data stored in the 3D structure 
	  //for the ray tracer by means of pointers.
	  //The structure acts as a  hash table - the remainined (%) of dividing
	  //the cluster voxel ID against NUMBEROFCLUSTERVOXELMAPS indicates which map
	  //to look up. This is chosen because the cluster voxel IDs tend to be spare
	  //but large blocks of consecutive numbers.
	  std::vector<
	    std::map<site_t, SiteData_t*> > 
	    mClusterVoxelDataPointers;


	  const geometry::LatticeData*& mLatticeData;
	  short int *mClusterIdOfBlock;
      
	  static const unsigned int NUMBEROFCLUSTERVOXELMAPS = 1024;

	  static const short int NOTASSIGNEDTOCLUSTER = -1; 

	  static const Location<site_t> mNeighbours[26];
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


	void UpdateRayData(const SiteData_t* iSiteData,
			   float ray_t,
			   float ray_segment,
			   Ray* bCurrentRay);

	void TraverseVoxels(const float block_min[3],
			    const float block_x[3],
			    const SiteData_t* iSiteData,
			    float t,
			    Ray* bCurrentRay,
			    const bool xyz_is_1[3]);

	void TraverseBlocks(const Cluster* cluster,
			    const bool xyz_Is_1[3],
			    const float ray_dx[3],
			    Ray *bCurrentRay);

	void AABBvsRay(const AABB* aabb,
		       const float inverseDirection[3],
		       const bool xyzComponentIsPositive[3],
		       float* t_near,
		       float* t_far);

	void UpdateColour(float dt, const float palette[3], float col[3]);

	ClusterBuilder mClusterBuilder;
	
	const geometry::LatticeData* mLatDat;

	const DomainStats* mDomainStats;
	Screen* mScreen;
	Viewpoint* mViewpoint;
	VisSettings* mVisSettings;

	const float mBlockSizeFloat;
	const float mBlockSizeInverse;
	const site_t block_size2, block_size3, block_size_1;
	const site_t blocksYz;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_H
