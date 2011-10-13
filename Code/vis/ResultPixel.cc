#include "vis/ResultPixel.h"

namespace hemelb
{
  namespace vis
  {
    ResultPixel::ResultPixel(const BasicPixel* glyph) :
      BasicPixel(glyph->GetI(), glyph->GetJ()), hasGlyph(true), rayPixel(NULL), streakPixel(NULL)
    {

    }

    ResultPixel::ResultPixel(const raytracer::RayPixel* ray) :
      BasicPixel(ray->GetI(), ray->GetJ()), hasGlyph(false), rayPixel(ray), streakPixel(NULL)
    {

    }

    ResultPixel::ResultPixel(const StreakPixel* streak) :
      BasicPixel(streak->GetI(), streak->GetJ()), hasGlyph(false), rayPixel(NULL),
          streakPixel(streak)
    {

    }

    const raytracer::RayPixel* ResultPixel::GetRayPixel() const
    {
      return rayPixel;
    }

    void ResultPixel::Combine(const ResultPixel& other)
    {
      if (other.hasGlyph)
      {
        hasGlyph = true;
      }

      if (other.rayPixel != NULL)
      {
        rayPixel = other.rayPixel;
      }

      if (other.streakPixel != NULL)
      {
        streakPixel = other.streakPixel;
      }
    }

    void ResultPixel::WritePixel(int *pixel_index,
                                 unsigned char rgb_data[12],
                                 const DomainStats* iDomainStats,
                                 const VisSettings* visSettings) const
    {
      const int bits_per_char = sizeof(char) * 8;
      *pixel_index = (i << (2 * bits_per_char)) + j;

      if (rayPixel != NULL)
      {
        // store velocity volume rendering colour
        float dt = rayPixel->GetDT();
        const float* velArray = rayPixel->GetVelArray();
        MakePixelColour(int (255.0F * velArray[0] / dt),
                        int (255.0F * velArray[1] / dt),
                        int (255.0F * velArray[2] / dt),
                        &rgb_data[0]);

        const float* stressArr = rayPixel->GetStressArray();
        float stress = rayPixel->GetStress();

        if (visSettings->mStressType != lb::ShearStress)
        {
          // store von Mises stress volume rendering colour
          MakePixelColour(int (stressArr[0] / dt),
                          int (stressArr[1] / dt),
                          int (stressArr[2] / dt),
                          &rgb_data[3]);
        }
        else if (stress < ((float) NO_VALUE))
        {
          float stress_col[3];
          PickColour(stress, stress_col);

          // store wall shear stress colour
          MakePixelColour(int (255.0F * stress_col[0]), int (255.0F * stress_col[1]), int (255.0F
              * stress_col[2]), &rgb_data[3]);
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

      float density = rayPixel == NULL
        ? 0.0F
        : rayPixel->GetDensity();
      float stress = rayPixel == NULL
        ? 0.0F
        : rayPixel->GetStress();

      if (visSettings->mStressType != lb::ShearStress && visSettings->mode
          == VisSettings::ISOSURFACES)
      {
        float density_col[3], stress_col[3];
        PickColour(density, density_col);
        PickColour(stress, stress_col);

        // store wall pressure colour
        MakePixelColour(int (255.0F * density_col[0]), int (255.0F * density_col[1]), int (255.0F
            * density_col[2]), &rgb_data[6]);

        // store von Mises stress colour
        MakePixelColour(int (255.0F * stress_col[0]), int (255.0F * stress_col[1]), int (255.0F
            * stress_col[2]), &rgb_data[9]);

      }
      else if (visSettings->mStressType != lb::ShearStress && visSettings->mode
          == VisSettings::ISOSURFACESANDGLYPHS)
      {
        float density_col[3], stress_col[3];
        PickColour(density, density_col);
        PickColour(stress, stress_col);

        if (rayPixel != NULL)
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
          MakePixelColour(int (127.5F * density_col[0]), int (127.5F * density_col[1]), int (127.5F
              * density_col[2]), &rgb_data[6]);

          // store von Mises stress (+glyph) colour
          MakePixelColour(int (127.5F * stress_col[0]), int (127.5F * stress_col[1]), int (127.5F
              * stress_col[2]), &rgb_data[9]);
        }
        else
        {
          for (int ii = 6; ii < 12; ++ii)
          {
            rgb_data[ii] = 0;
          }
        }

      }
      else if (streakPixel != NULL)
      {
        float scaled_vel = (float) (streakPixel->GetParticleVelocity()
            * iDomainStats->velocity_threshold_max_inv);
        float particle_col[3];
        PickColour(scaled_vel, particle_col);

        // store particle colour
        MakePixelColour(int (255.0F * particle_col[0]), int (255.0F * particle_col[1]), int (255.0F
            * particle_col[2]), &rgb_data[6]);

        for (int ii = 9; ii < 12; ++ii)
        {
          rgb_data[ii] = rgb_data[ii - 3];
        }
      }

      else
      {
        // store pressure colour
        rgb_data[6] = rgb_data[7] = rgb_data[8]
            = (unsigned char) util::NumericalFunctions::enforceBounds(int (127.5F * density),
                                                                      0,
                                                                      127);

        // store shear stress or von Mises stress
        if (stress < ((float) NO_VALUE))
        {
          rgb_data[9] = rgb_data[10] = rgb_data[11]
              = (unsigned char) util::NumericalFunctions::enforceBounds(int (127.5F * stress),
                                                                        0,
                                                                        127);
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
      colour[1] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * (float) fabs(value
          - 0.5F), 0.F, 1.F);
      colour[2] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * value, 0.F, 1.F);
    }

    void ResultPixel::MakePixelColour(int rawRed, int rawGreen, int rawBlue, unsigned char* dest)
    {
      dest[0] = (unsigned char) util::NumericalFunctions::enforceBounds(rawRed, 0, 255);
      dest[1] = (unsigned char) util::NumericalFunctions::enforceBounds(rawGreen, 0, 255);
      dest[2] = (unsigned char) util::NumericalFunctions::enforceBounds(rawBlue, 0, 255);
    }
  }
}
