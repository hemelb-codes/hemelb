#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis 
  {
    namespace raytracer 
    {
      RayTracer::ClusterBuilder::VolumeTraverser::VolumeTraverser()
	: mCurrentLocation(0)
      {
	mCurrentNumber=0;
      }
	
      Location<site_t> RayTracer::ClusterBuilder::VolumeTraverser::GetCurrentLocation()
      {
	return mCurrentLocation;
      }

      site_t  RayTracer::ClusterBuilder::VolumeTraverser::GetCurrentIndex()
      {
	return mCurrentNumber;
      }

      site_t RayTracer::ClusterBuilder::VolumeTraverser::
      GetIndexFromLocation(Location<site_t> iLocation)
      {
	return ((iLocation.x * GetYCount() + iLocation.y) 
		* GetZCount()) + iLocation.z;
      }
	    

      bool RayTracer::ClusterBuilder::VolumeTraverser::TraverseOne()
      {
	mCurrentNumber++;

	mCurrentLocation.z++;
	if(mCurrentLocation.z < GetZCount())
	{
	  return true;
	}
	
		
	mCurrentLocation.z = 0;
	mCurrentLocation.y++;
	if(mCurrentLocation.y < GetYCount())
	{
	  return true;
	}
		
	mCurrentLocation.y = 0;
	mCurrentLocation.x++;
			 
	if(mCurrentLocation.x < GetXCount())
	{
	  return true;
	}
	return false;
      }
    }
  }
}
