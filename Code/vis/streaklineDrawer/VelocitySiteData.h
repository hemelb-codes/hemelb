// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H
#define HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H

#include "util/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      // Class to hold information about the velocity field at some point.
      class VelocitySiteData
      {
        public:
          VelocitySiteData() :
              counter(-1), proc_id(-1), site_id(-1), velocity(NO_VALUE)
          {
          }

          /**
           * Counter describes how many iterations have passed since the objects creation.
           * It allows for quickly checking whether the local velocity field has been
           * calculated this iteration.
           **/
          site_t counter;

          /**
           * The rank this volume of the geometry lives on.
           */
          proc_t proc_id;

          /**
           * The site id of the volume represented.
           */
          site_t site_id;

          /**
           * The velocity calculated for that site.
           */
          util::Vector3D<float> velocity;
      };
    }
  }
}

#endif //HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H
