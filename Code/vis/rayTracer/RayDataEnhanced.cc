// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "vis/rayTracer/RayDataEnhanced.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      namespace DepthCuing
      {
        const float Mist::SurfaceNormalLightnessRange = 0.3F;
        const float Mist::ParallelSurfaceAttenuation = 0.5F;
        const float Mist::LowestLightness = 0.3F;
        const float Mist::VelocityHueMin = 240.0F;
        const float Mist::VelocityHueRange = 120.0F;
        const float Mist::VelocitySaturation = 1.0F;
        const float Mist::StressHue = 230.0F;
        const float Mist::StressSaturationRange = 0.5F;
        const float Mist::StressSaturationMin = 0.5F;

        const float None::SurfaceNormalLightnessRange = 0.5F;
        const float None::ParallelSurfaceAttenuation = 0.75F;
        const float None::LowestLightness = 0.3F;
        const float None::VelocityHueMin = 240.0F;
        const float None::VelocityHueRange = 120.0F;
        const float None::VelocitySaturation = 1.0F;
        const float None::StressHue = 230.0F;
        const float None::StressSaturationRange = 0.5F;
        const float None::StressSaturationMin = 0.5F;

        const float Darkness::SurfaceNormalLightnessRange = 0.3F;
        const float Darkness::ParallelSurfaceAttenuation = 0.5F;
        const float Darkness::LowestLightness = 0.0F;
        const float Darkness::VelocityHueMin = 240.0F;
        const float Darkness::VelocityHueRange = 120.0F;
        const float Darkness::VelocitySaturation = 1.0F;
        const float Darkness::StressHue = 230.0F;
        const float Darkness::StressSaturationRange = 0.5F;
        const float Darkness::StressSaturationMin = 0.5F;
      }
    }
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced>::RegisterMpiDataType()
  {
    MPI_Datatype ret = vis::raytracer::RayDataEnhanced::GetMpiType();
    MPI_Type_commit(&ret);
    return ret;
  }
}
