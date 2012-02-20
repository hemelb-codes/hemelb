#ifndef HEMELB_VIS_RAYTRACER_SITEDATA_H
#define HEMELB_VIS_RAYTRACER_SITEDATA_H

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
          float density;
          float velocity;
          float stress;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_SITEDATA_H
