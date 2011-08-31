#ifndef HEMELB_VIS_RAYTRACER_CLUSTERSHARED_H
#define HEMELB_VIS_RAYTRACER_CLUSTERSHARED_H

#include <vector>

#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/SiteData.h"
#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
	//The cluster structure stores data relating to the clusters
	//used by the RayTracer, in an optimal format
	//Cluster are produced by the ClusterSharedFactory
	//Caution: the data within the flow field is altered by means
	//of pointers obtained from the GetClusterSharedVoxelDataPointer
	//method
      template <typename T>
	class ClusterShared : public Cluster<T>
	{
	public:
	  ClusterShared() {}

	  unsigned int GetBlockIdFrom3DBlockLocation(
	    Vector3D<unsigned int>iLocation) const
	  {
	    return iLocation.x*blocksY*blocksZ
	      + iLocation.y*blocksZ
	      + iLocation.z;
	  }

	  //Returns true if there is site data for a given block
	  bool BlockContainsSites(site_t iBlockNumber) const
	  {
	    return !SiteData[iBlockNumber].empty();
	  }

	  void ResizeVectorsForBlockShared(site_t iBlockNumber, site_t iSize)
	  {
	    //By default all values are -1 for solids
	    SiteData[iBlockNumber].resize(iSize, SiteData_t(-1.0F)); 
	  }
	  
	  //Get SiteData arary for site
	  const SiteData_t* GetSiteData(site_t iBlockNumber) const
	  {
	    return &(SiteData[iBlockNumber][0]);
	  }	    

	  const SiteData_t* GetSiteData(site_t iBlockNumber, site_t iSiteNumber) const
	  {
	    return &(SiteData[iBlockNumber][iSiteNumber]);
	  }
      
	  static bool NeedsWallNormals() //overridden
	  {
	    return false; // By default
	  }

	  //The min and maximum site location, in site units 
	  //relative to the centre of the lattice
	  Vector3D<float> minSite;
	  Vector3D<float> maxSite;

	  //Stores the lowest x, y and z block location of the ClusterShared 
	  //in terms of site units relative to the centre location
	  Vector3D<float> minBlock;
        
	  //Stores the size of the cluster in terms of the number of blocks
	  unsigned short int blocksX;
	  unsigned short int blocksY;
	  unsigned short int blocksZ;

	  std::vector<std::vector<SiteData_t> > SiteData;

	protected:
	  void ResizeSharedVectors()
	  {
	    SiteData.resize(blocksX*blocksY*blocksZ);
	  }
	  
	};

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERSHARED_H
