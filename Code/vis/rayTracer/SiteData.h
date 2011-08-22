#ifndef HEMELB_VIS_SITEDATA_H
#define HEMELB_VIS_SITEDATA_H

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      //Stores the data about an individual voxel 
      struct SiteData_t
      {
      public:
	float Density;
	float Velocity;
	float Stress;

      SiteData_t(float iValue) :
	Density(iValue), Velocity(iValue), Stress(iValue) {}
      };
    }
  }
}

#endif // HEMELB_VIS_SITEDATA_H
