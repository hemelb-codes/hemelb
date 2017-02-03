// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "vis/ResultPixel.h"

namespace hemelb
{
  namespace vis
  {
    ResultPixel::ResultPixel(const BasicPixel* glyph) :
        BasicPixel(glyph->GetI(), glyph->GetJ()), hasGlyph(true), normalRayPixel(nullptr),
            streakPixel(nullptr)
    {

    }

    ResultPixel::ResultPixel(const raytracer::RayDataNormal* ray) :
        BasicPixel(ray->GetI(), ray->GetJ()), hasGlyph(false), normalRayPixel(ray),
            streakPixel(nullptr)
    {

    }

    ResultPixel::ResultPixel(const streaklinedrawer::StreakPixel* streak) :
        BasicPixel(streak->GetI(), streak->GetJ()), hasGlyph(false), normalRayPixel(nullptr),
            streakPixel(streak)
    {

    }

    const raytracer::RayDataNormal* ResultPixel::GetRayPixel() const
    {
      return normalRayPixel;
    }

    void ResultPixel::Combine(const ResultPixel& other)
    {
      if (other.hasGlyph)
      {
        hasGlyph = true;
      }

      if (other.normalRayPixel != nullptr)
      {
        normalRayPixel = other.normalRayPixel;
      }

      if (other.streakPixel != nullptr)
      {
        streakPixel = other.streakPixel;
      }
    }

    void ResultPixel::WritePixel(unsigned* pixel_index, unsigned char rgb_data[12],
                                 const DomainStats& iDomainStats,
                                 const VisSettings& visSettings) const
    {
      const int bits_per_char = sizeof(char) * 8;
      *pixel_index = (i << (2 * bits_per_char)) + j;

      if (normalRayPixel != nullptr)
      {
        // store velocity volume rendering colour
        normalRayPixel->GetVelocityColour(rgb_data, visSettings, iDomainStats);

        float stress = normalRayPixel->GetNearestStress();

        if (visSettings.mStressType != lb::ShearStress)
        {
          normalRayPixel->GetStressColour(&rgb_data[3], visSettings, iDomainStats);
        }
        else if (stress < (float) NO_VALUE)
        {
          float stress_col[3];
          PickColour(stress, stress_col);

          // store wall shear stress colour
          MakePixelColour(int(255.0F * stress_col[0]),
                          int(255.0F * stress_col[1]),
                          int(255.0F * stress_col[2]),
                          &rgb_data[3]);
        }
        else
        {
          rgb_data[3] = rgb_data[4] = rgb_data[5] = 0;
        }
      }
      else
      {
        for (int ii = 0; ii < 6; ++ii)
        {
          rgb_data[ii] = 255;
        }
      }

      float density = normalRayPixel == nullptr ?
        0.0F :
        normalRayPixel->GetNearestDensity();
      float stress = normalRayPixel == nullptr ?
        0.0F :
        normalRayPixel->GetNearestStress();

      if (visSettings.mStressType != lb::ShearStress
          && visSettings.mode == VisSettings::ISOSURFACES)
      {
        float density_col[3], stress_col[3];
        PickColour(density, density_col);
        PickColour(stress, stress_col);

        // store wall pressure colour
        MakePixelColour(int(255.0F * density_col[0]),
                        int(255.0F * density_col[1]),
                        int(255.0F * density_col[2]),
                        &rgb_data[6]);

        // store von Mises stress colour
        MakePixelColour(int(255.0F * stress_col[0]),
                        int(255.0F * stress_col[1]),
                        int(255.0F * stress_col[2]),
                        &rgb_data[9]);

      }
      else if (visSettings.mStressType != lb::ShearStress
          && visSettings.mode == VisSettings::ISOSURFACESANDGLYPHS)
      {
        float density_col[3], stress_col[3];
        PickColour(density, density_col);
        PickColour(stress, stress_col);

        if (normalRayPixel != nullptr)
        {
          if (!hasGlyph)
          {
            density_col[0] += 1.0F;
            density_col[1] += 1.0F;
            density_col[2] += 1.0F;

            stress_col[0] += 1.0F;
            stress_col[1] += 1.0F;
            stress_col[2] += 1.0F;
          }

          // store wall pressure (+glyph) colour
          MakePixelColour(int(127.5F * density_col[0]),
                          int(127.5F * density_col[1]),
                          int(127.5F * density_col[2]),
                          &rgb_data[6]);

          // store von Mises stress (+glyph) colour
          MakePixelColour(int(127.5F * stress_col[0]),
                          int(127.5F * stress_col[1]),
                          int(127.5F * stress_col[2]),
                          &rgb_data[9]);
        }
        else
        {
          for (int ii = 6; ii < 12; ++ii)
          {
            rgb_data[ii] = 0;
          }
        }

      }
      else if (streakPixel != nullptr)
      {
        float scaled_vel = (float) (streakPixel->GetParticleVelocity()
            * iDomainStats.velocity_threshold_max_inv);
        float particle_col[3];
        PickColour(scaled_vel, particle_col);

        // store particle colour
        MakePixelColour(int(255.0F * particle_col[0]),
                        int(255.0F * particle_col[1]),
                        int(255.0F * particle_col[2]),
                        &rgb_data[6]);

        for (int ii = 9; ii < 12; ++ii)
        {
          rgb_data[ii] = rgb_data[ii - 3];
        }
      }
      else
      {
        // store pressure colour
        rgb_data[6] = rgb_data[7] = rgb_data[8] =
            (unsigned char) util::NumericalFunctions::enforceBounds(int(127.5F * density), 0, 127);

        // store shear stress or von Mises stress
        if (stress < ((float) NO_VALUE))
        {
          rgb_data[9] = rgb_data[10] = rgb_data[11] =
              (unsigned char) util::NumericalFunctions::enforceBounds(int(127.5F * stress), 0, 127);
        }
        else
        {
          rgb_data[9] = rgb_data[10] = rgb_data[11] = 0;
        }
      }
    }

    void ResultPixel::PickColour(float value, float colour[3])
    {
      colour[0] = util::NumericalFunctions::enforceBounds<float>(4.F * value - 2.F, 0.F, 1.F);
      colour[1] = util::NumericalFunctions::enforceBounds<float>(2.F
                                                                     - 4.F
                                                                         * (float) fabs(value
                                                                             - 0.5F),
                                                                 0.F,
                                                                 1.F);
      colour[2] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * value, 0.F, 1.F);
    }

    void ResultPixel::MakePixelColour(int rawRed, int rawGreen, int rawBlue, unsigned char* dest)
    {
      dest[0] = (unsigned char) util::NumericalFunctions::enforceBounds(rawRed, 0, 255);
      dest[1] = (unsigned char) util::NumericalFunctions::enforceBounds(rawGreen, 0, 255);
      dest[2] = (unsigned char) util::NumericalFunctions::enforceBounds(rawBlue, 0, 255);
    }

    void ResultPixel::LogDebuggingInformation() const
    {
      log::Logger::Log<log::Trace, log::OnePerCore>("Pixel at (%i,%i) with (ray,streak,glyph)=(%i,%i,%i)",
                                                    GetI(),
                                                    GetJ(),
                                                    normalRayPixel != nullptr,
                                                    streakPixel != nullptr,
                                                    hasGlyph);

      if (normalRayPixel != nullptr)
      {
        normalRayPixel->LogDebuggingInformation();
      }

      if (streakPixel != nullptr)
      {
        streakPixel->LogDebuggingInformation();
      }
    }

  }
}
