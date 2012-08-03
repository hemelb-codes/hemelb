// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
