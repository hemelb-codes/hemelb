// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_STEERING_COMMON_TAGS_H
#define HEMELB_STEERING_COMMON_TAGS_H

#include "steering/common/Tag.h"
#include "steering/common/Steerer.h"

namespace hemelb
{
  namespace steering
  {
    HEMELB_STEERING_DECLARE_TAG(ctr_x, float);
    HEMELB_STEERING_DECLARE_TAG(ctr_y, float);
    HEMELB_STEERING_DECLARE_TAG(ctr_z, float);

    HEMELB_STEERING_DECLARE_TAG(longitude, float);
    HEMELB_STEERING_DECLARE_TAG(latitude, float);

    HEMELB_STEERING_DECLARE_TAG(zoom, float);
    HEMELB_STEERING_DECLARE_TAG(brightness, float);

    HEMELB_STEERING_DECLARE_TAG(velocity_max, float);
    HEMELB_STEERING_DECLARE_TAG(stress_max, float);

    HEMELB_STEERING_DECLARE_TAG(pressure_min, float);
    HEMELB_STEERING_DECLARE_TAG(pressure_max, float);

    HEMELB_STEERING_DECLARE_TAG(glyph_length, float);

    HEMELB_STEERING_DECLARE_TAG(pixels_x, unsigned int);
    HEMELB_STEERING_DECLARE_TAG(pixels_y, unsigned int);

    HEMELB_STEERING_DECLARE_TAG(terminate_signal, bool);
    HEMELB_STEERING_DECLARE_TAG(vis_mode, int);
    HEMELB_STEERING_DECLARE_TAG(streaklines_per_pulsatile_period, unsigned int);
    HEMELB_STEERING_DECLARE_TAG(streakline_length, float);
    HEMELB_STEERING_DECLARE_TAG(doRendering, bool);

  }
}

#endif // HEMELB_STEERING_COMMON_TAGS_H
