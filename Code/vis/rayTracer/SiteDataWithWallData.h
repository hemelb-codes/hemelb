#ifndef HEMELB_VIS_SITEDATAWITHWALLDATA_H
#define HEMELB_VIS_SITEDATAWITHWALLDATA_H

#include "vis/rayTracer/SiteData.h";

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      //Stores the data about an individual voxel 
      struct SiteDataWithWallData : public SiteData
      {
      public:
	double* WallData;

      SiteDataWithWallData(float iValue) :
	SiteData(iValue)
	{
	  WallData = NULL;
	}
      };
    }
  }
}

#endif // HEMELB_VIS_SITEDATAWITHWALLDATA_H
