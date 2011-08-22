#ifndef HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H
#define HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H

#include <stack>

#include "vis/Vector3D.h"
#include "vis/rayTracer/BlockTraverser.h"
#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class ClusterBuilder
      {
      public:
	ClusterBuilder
	  (const geometry::LatticeData*& iLatDat);
        ~ClusterBuilder();
	  
	void BuildClusters();
	  
	std::vector<Cluster*>& GetClusters();

	SiteData_t* GetClusterVoxelDataPointer(site_t iSiteId);

      protected:
	void UpdateDensityVelocityAndStress(site_t iBlockNum, unsigned int iClusterId, unsigned int iSiteIdOnBlock, unsigned int iClusterVoxelSiteId);

	//Caution: the data within mClusters is altered by means
	//of pointers obtained from the GetClusterVoxelDataPointer
	//method. No insertion of copying must therefore take place
	//on mClusters once building is complete
	std::vector<Cluster*> mClusters;

	const geometry::LatticeData*& mLatticeData;

      private:
	//Locates all the clusters in the lattice structure and the
	void LocateClusters();
 
	//Locates all the clusters in the lattice structure and the
	void FindNewCluster();

	//Adds neighbouring blocks of the input location to the input stack
	void AddNeighbouringBlocks(Vector3D<site_t> iCurrentLocation,
				   std::stack<Vector3D<site_t> >& ioBlocksToVisit);

	//Returns true if there are sites in the given block associated with the
	//local processor rank
	bool AreSitesAssignedToLocalProcessorRankInBlock
	  (geometry::LatticeData::BlockData * iBlock);

	//Creates a cluster of the correct type
	virtual Cluster* CreateNewCluster();

	//Adds a new cluster by taking in the required data in interget format
	//and converting it to that used by the raytracer
	//NB: Futher processing is required on the cluster before it can be used
	//by the ray tracer, which is handled by the ProcessCluster method
	void AddCluster(Vector3D<site_t> iClusterBlockMin,
			Vector3D<site_t> iClusterBlockMax,
			Vector3D<site_t> iClusterVoxelMin,
			Vector3D<site_t> iClusterVoxelMax);

	//Adds "flow-field" data to the cluster
	void ProcessCluster(unsigned int iClusterId);

	void UpdateSiteData
	  (site_t iBlockId, site_t iBlockNum,  unsigned int iClusterId,
	   Vector3D<site_t>i_block_coordinates);

	virtual void UpdateSiteDataAtSite
	  (site_t iBlockId, site_t iBlockNum, 
	   unsigned int iClusterId, unsigned int iSiteIdOnBlock);
	
	Vector3D<site_t> GetSiteCoordinatesOfBlock
	  (site_t iClusterId, Vector3D<site_t> offset);

	SiteData_t* GetDataPointerClusterVoxelSiteId(site_t iSiteId);

	void SetDataPointerForClusterVoxelSiteId
	  (site_t iClusterVortexSiteId,
	   SiteData_t* iDataPointer);

	BlockTraverser mBlockTraverser;

	std::vector<Vector3D<site_t> > mClusterBlockMins;

	//This allows a cluster voxel site ID (as part of the 1D structure for)
	//storing sites to be mapped to the data stored in the 3D structure
	//for the ray tracer by means of pointers.
	std::vector<SiteData_t*> mClusterVoxelDataPointers;

	short int *mClusterIdOfBlock;

	static const short int NOTASSIGNEDTOCLUSTER = -1;

	static const Vector3D<site_t> mNeighbours[26];
      };
      
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H
