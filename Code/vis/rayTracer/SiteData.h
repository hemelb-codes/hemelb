// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
