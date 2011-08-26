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

      void VolumeTraverser::SetCurrentLocation(
	const Vector3D<site_t>& iLocation)
      {
	mCurrentLocation = iLocation;
	mCurrentNumber = GetIndexFromLocation(iLocation);
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

      void VolumeTraverser::IncrementX()
      {
	mCurrentLocation.x++;
	mCurrentNumber += GetZCount() * GetYCount();
      }

      void VolumeTraverser::IncrementY()
      {
	mCurrentLocation.y++;
	mCurrentNumber += GetZCount();
      }

       void VolumeTraverser::IncrementZ()
      {
	mCurrentLocation.z++;
	mCurrentNumber++;
      }

      void VolumeTraverser::DecrementX()
      {
	mCurrentLocation.x--;
	mCurrentNumber -= GetZCount() * GetYCount();
      }

      void VolumeTraverser::DecrementY()
      {
	mCurrentLocation.y--;
	mCurrentNumber -= GetZCount();
      }

       void VolumeTraverser::DecrementZ()
      {
	mCurrentLocation.z--;
	mCurrentNumber--;
      }
      
      bool VolumeTraverser::CurrentLocationValid()
      {
	if (GetCurrentIndex() < 0)
	{
	  return false;
	}

	if (mCurrentLocation.x < 0 ||
	    mCurrentLocation.y < 0 ||
	    mCurrentLocation.z < 0 ||
	    mCurrentLocation.x >= GetXCount() ||
	    mCurrentLocation.y >= GetYCount() ||
	    mCurrentLocation.z >= GetZCount())
	{
	  return false;
	}
	
	return true;
      }
    }
  }
}
