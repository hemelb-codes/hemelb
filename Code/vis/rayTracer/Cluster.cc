#include "vis/rayTracer/Cluster.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      Cluster::Cluster() { }

      unsigned int Cluster::GetBlockIdFrom3DBlockLocation(
	Vector3D<unsigned int>iLocation) const
      {
	return iLocation.x*blocksY*blocksZ
	  + iLocation.y*blocksZ
	  + iLocation.z;
      }

      void Cluster::ResizeVectors()
      {
	SiteData.resize(blocksX*blocksY*blocksZ);
      }

      bool Cluster::BlockContainsSites(site_t  iBlockNumber) const
      {
	return !SiteData[iBlockNumber].empty();
      }
      
      const SiteData_t* Cluster::GetSiteData(site_t iBlockNumber) const
      {
	return &(SiteData[iBlockNumber][0]);
      }

      const SiteData_t* Cluster::GetSiteData(site_t iBlockNumber, site_t iSiteNumber) const
      {
	return &(SiteData[iBlockNumber][iSiteNumber]);
      }
    }
  }
}
