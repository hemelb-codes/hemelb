#ifndef HEMELB_VIS_SITEDATAWITHWALLDATA_H
#define HEMELB_VIS_SITEDATAWITHWALLDATA_H

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      //Stores the data about an individual voxel 
      struct SiteDataWithWallData : public SiteDate
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
