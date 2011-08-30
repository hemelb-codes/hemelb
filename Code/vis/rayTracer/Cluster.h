#ifndef HEMELB_VIS_RAYTRACER_CLUSTER_H
#define HEMELB_VIS_RAYTRACER_CLUSTER_H

#include <vector>

#include "vis/rayTracer/SiteData.h"
#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
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

	  unsigned int GetBlockIdFrom3DBlockLocation(
	    Vector3D<unsigned int>iLocation) const;

	  //Resizes the vectors so as to be the correct size based on the stored sizes
	  virtual void ResizeVectors();

	  //Returns true if there is site data for a given block
	  bool BlockContainsSites(site_t iBlockNumber) const;
	  
	  //Get SiteData arary for site
	  const SiteData_t* GetSiteData(site_t iBlockNumber) const;
	  
	  //The min and maximum site location, in site units 
	  //relative to the centre of the lattice
	  Vector3D<float> minSite;
	  Vector3D<float> maxSite;

	  //Stores the lowest x, y and z block location of the Cluster 
	  //in terms of site units relative to the centre location
	  Vector3D<float> minBlock;
        
	  //Stores the size of the cluster in terms of the number of blocks
	  unsigned short int blocksX;
	  unsigned short int blocksY;
	  unsigned short int blocksZ;

	  std::vector<std::vector<SiteData_t> > SiteData;

	};

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTER_H
