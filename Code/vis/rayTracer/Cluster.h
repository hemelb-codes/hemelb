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
      //used by the RayTracer, in an optimal format
      //Cluster are produced by the ClusterSharedFactory
      //Caution: the data within the flow field is altered by means
      //of pointers obtained from the GetClusterSharedVoxelDataPointer
      //method
      template <typename Derived>
	class Cluster
      {
      public:
	unsigned int GetBlockIdFrom3DBlockLocation(
	  Vector3D<unsigned int>iLocation) const
	{
	  return static_cast<Derived*>(this)->GetBlockIdFrom3DBlockLocation(iLocation);
	}

	//Resizes the vectors so as to be the correct size based on the stored sizes
	void ResizeVectors()
	{
	  static_cast<Derived*>(this)->ResizeVectors();
	}

	void ResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
	{
	  static_cast<Derived*>(this)->ResizeVectorsForBlock(iBlockNumber, iSize);
	}

	//Returns true if there is site data for a given block
	bool BlockContainsSites(site_t iBlockNumber) const
	{
	  return static_cast<Derived*>(this)->BlockContainsSites(iBlockNumber);
	}
	  
	//Get SiteData arary for site
	const SiteData_t* GetSiteData(site_t iBlockNumber) const
	{
	  return static_cast<Derived*>(this)->GetSiteData(iBlockNumber);
	}	    

	const SiteData_t* GetSiteData(site_t iBlockNumber, site_t iSiteNumber) const
	{
	  return static_cast<Derived*>(this)->GetSiteData(iBlockNumber, iSiteNumber);
	}
      
	double const* GetWallData(site_t iBlockNumber, site_t iSiteNumber) const
	{
	  return static_cast<Derived*>(this)->GetWallData(iBlockNumber, iSiteNumber);
	}	  

	void SetWallData(site_t iBlockNumber, site_t iSiteNumber, double* iData)
	{
	  return static_cast<Derived*>(this)->SetWallData(iBlockNumber, iSiteNumber, iData);
	}

	static bool NeedsWallNormals()
	{
	  return Derived::NeedsWallNormals();
	}
      };

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTER_H
