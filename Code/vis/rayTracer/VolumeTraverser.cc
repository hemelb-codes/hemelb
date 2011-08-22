#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis 
  {
    namespace raytracer 
    {
      VolumeTraverser::VolumeTraverser()
	: mCurrentLocation(0)
      {
	mCurrentNumber=0;
      }
	
      Vector3D<site_t> VolumeTraverser::GetCurrentLocation()
      {
	return mCurrentLocation;
      }

      site_t  VolumeTraverser::GetCurrentIndex()
      {
	return mCurrentNumber;
      }

      site_t VolumeTraverser::
      GetIndexFromLocation(Vector3D<site_t> iLocation)
      {
	return ((iLocation.x * GetYCount() + iLocation.y) 
		* GetZCount()) + iLocation.z;
      }
	    

      bool VolumeTraverser::TraverseOne()
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
