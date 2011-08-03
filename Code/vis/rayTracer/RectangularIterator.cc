#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
    namespace vis 
    {
	namespace raytracer 
	{
	    RayTracer::ClusterBuilder::RectangularIterator::RectangularIterator()
		: mCurrentLocation(0)
	    {
		mCurrentNumber=0;
	    }
	
	    Location RayTracer::ClusterBuilder::RectangularIterator::GetLocation()
	    {
		return mCurrentLocation;
	    }

	    site_t  RayTracer::ClusterBuilder::RectangularIterator::CurrentNumber()
	    {
		return mCurrentNumber;
	    }

	    site_t RayTracer::ClusterBuilder::RectangularIterator::
	    GetNumberFromLocation(Location iLocation)
	    {
		return ((iLocation.i * GetYCount() + iLocation.j) 
			* GetZCount()) + iLocation.k;
	    }
	    

	    bool RayTracer::ClusterBuilder::RectangularIterator::Iterate()
	    {
		mCurrentNumber++;

		mCurrentLocation.k++;
		if(mCurrentLocation.k < GetZCount())
		{
		    return true;
		}
	
		
		mCurrentLocation.k = 0;
		mCurrentLocation.j++;
		if(mCurrentLocation.j < GetYCount())
		{
		    return true;
		}
		
		mCurrentLocation.j = 0;
		mCurrentLocation.i++;
			 
		if(mCurrentLocation.i < GetXCount())
		{
		    return true;
		}
		return false;
	    }
 	}
    }
}
