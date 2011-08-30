#include <math.h>

#include "vis/ColPixel.h"
#include "util/utilityFunctions.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {

    ColPixel::PixelId::PixelId(unsigned int i_, unsigned int j_) :
      isRt(false), isGlyph(false), isStreakline(false), i(i_), j(j_)
    {
    }

    ColPixel::PixelId::PixelId() :
      isRt(false), isGlyph(false), isStreakline(false), i(0), j(0)
    {
    }

    ColPixel::ColPixel()
    {
    }

    ColPixel::ColPixel(int iIn, int jIn) :
      i(iIn, jIn)
    {
      i.isGlyph = true;
    }

    ColPixel::ColPixel(int iIn,
                       int jIn,
                       float particleVelocity,
                       float particleZ,
                       int particleInletId) :
      i(iIn, jIn)
    {
      i.isStreakline = true;

      particle_vel = particleVelocity;
      particle_z = particleZ;
      particle_inlet_id = particleInletId;
    }

    ColPixel::ColPixel(int iIn,
                       int jIn,
                       float tIn,
                       float dtIn,
                       float densityIn,
                       float stressIn,
                       const float velocityColour[3],
                       const float stressColour[3]) :
      i(iIn, jIn)
    {
      i.isRt = true;

      LengthBeforeRayFirstCluster = tIn;
      dt = dtIn;

      density = densityIn;
      stress = stressIn;

      vel_r = velocityColour[0] * 255.0F;
      vel_g = velocityColour[1] * 255.0F;
      vel_b = velocityColour[2] * 255.0F;

      stress_r = stressColour[0] * 255.0F;
      stress_g = stressColour[1] * 255.0F;
      stress_b = stressColour[2] * 255.0F;
    }

    void ColPixel::MakePixelColour(int rawRed, int rawGreen, int rawBlue, unsigned char* dest) const
    {
      dest[0] = (unsigned char) util::NumericalFunctions::enforceBounds(rawRed, 0, 255);
      dest[1] = (unsigned char) util::NumericalFunctions::enforceBounds(rawGreen, 0, 255);
      dest[2] = (unsigned char) util::NumericalFunctions::enforceBounds(rawBlue, 0, 255);
    }

    /**
     * Merge data from the ColPixel argument into this pixel.
     */
    void ColPixel::MergeIn(const ColPixel *fromPixel, const VisSettings* visSettings)
    {
      // Merge raytracing data

      if (fromPixel->i.isRt && i.isRt)
      {
        // Both are ray-tracing
        vel_r += fromPixel->vel_r;
        vel_g += fromPixel->vel_g;
        vel_b += fromPixel->vel_b;

        if (visSettings->mStressType != lb::ShearStress)
        {
          stress_r += fromPixel->stress_r;
          stress_g += fromPixel->stress_g;
          stress_b += fromPixel->stress_b;
        }

        dt += fromPixel->dt;

        if (fromPixel->LengthBeforeRayFirstCluster < LengthBeforeRayFirstCluster)
        {
          LengthBeforeRayFirstCluster = fromPixel->LengthBeforeRayFirstCluster;
          density = fromPixel->density;
          stress = fromPixel->stress;
        }

      }
      else if (fromPixel->i.isRt && !i.isRt)
      {
        // Only the 'from' merge-pixel is ray-tracing
        vel_r = fromPixel->vel_r;
        vel_g = fromPixel->vel_g;
        vel_b = fromPixel->vel_b;

        if (visSettings->mStressType != lb::ShearStress)
        {
          stress_r = fromPixel->stress_r;
          stress_g = fromPixel->stress_g;
          stress_b = fromPixel->stress_b;
        }

        LengthBeforeRayFirstCluster = fromPixel->LengthBeforeRayFirstCluster;
        dt = fromPixel->dt;
        density = fromPixel->density;
        stress = fromPixel->stress;

        i.isRt = true;
      }
      // Done merging ray-tracing - (last combinations would be if from-pixel has no ray-tracing data)

      // Now merge glyph data
      if (visSettings->mStressType != lb::ShearStress && (visSettings->mode
          == VisSettings::ISOSURFACES || visSettings->mode == VisSettings::ISOSURFACESANDGLYPHS))
      {
        if (fromPixel->i.isGlyph)
        {
          i.isGlyph = true;
        }
      }
      else
      {
#ifndef NO_STREAKLINES
        // merge streakline data
        if (fromPixel->i.isStreakline)
        {
          if (!i.isStreakline)
          {
            particle_z = fromPixel->particle_z;
            particle_vel = fromPixel->particle_vel;
            particle_inlet_id = fromPixel->particle_inlet_id;

            i.isStreakline = true;
          }
          else
          {
            if (fromPixel->particle_z < particle_z)
            {
              particle_z = fromPixel->particle_z;
              particle_vel = fromPixel->particle_vel;
              particle_inlet_id = fromPixel->particle_inlet_id;
            }
          }
        }
#endif
      }
    }

    void ColPixel::rawWritePixel(int *pixel_index,
                                 unsigned char rgb_data[12],
                                 const DomainStats* iDomainStats,
                                 const VisSettings* visSettings) const
    {
      const int bits_per_char = sizeof(char) * 8;
      *pixel_index = (i.i << (2 * bits_per_char)) + i.j;

      if (i.isRt)
      {
        // store velocity volume rendering colour
        MakePixelColour(int(vel_r / dt), int(vel_g / dt), int(vel_b / dt), &rgb_data[0]);

        if (visSettings->mStressType != lb::ShearStress)
        {
          // store von Mises stress volume rendering colour
          MakePixelColour(int(stress_r / dt), int(stress_g / dt), int(stress_b / dt), &rgb_data[3]);
        }
        else if (stress < ((float) NO_VALUE))
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
      } // if (isRt)
      else
      {
        for (int ii = 0; ii < 6; ++ii)
        {
          rgb_data[ii] = 255;
        }
      }

      if (visSettings->mStressType != lb::ShearStress && visSettings->mode
          == VisSettings::ISOSURFACES)
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
      else if (visSettings->mStressType != lb::ShearStress && visSettings->mode
          == VisSettings::ISOSURFACESANDGLYPHS)
      {
        float density_col[3], stress_col[3];
        PickColour(density, density_col);
        PickColour(stress, stress_col);

        if (i.isRt)
        {
          if (!i.isGlyph)
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
      else if (i.isStreakline)
      {
        float scaled_vel = (float) (particle_vel * iDomainStats->velocity_threshold_max_inv);
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
        rgb_data[6] = rgb_data[7]
            = rgb_data[8] = (unsigned char) util::NumericalFunctions::enforceBounds(int(127.5F
                                                                                        * density),
                                                                                    0,
                                                                                    127);

        // store shear stress or von Mises stress
        if (stress < ((float) NO_VALUE))
        {
          rgb_data[9] = rgb_data[10] = rgb_data[11]
              = (unsigned char) util::NumericalFunctions::enforceBounds(int(127.5F * stress),
                                                                        0,
                                                                        127);
        }
        else
        {
          rgb_data[9] = rgb_data[10] = rgb_data[11] = 0;
        }
      }
    }

    void ColPixel::PickColour(float value, float colour[3])
    {
      colour[0] = util::NumericalFunctions::enforceBounds<float>(4.F * value - 2.F, 0.F, 1.F);
      colour[1] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * (float) fabs(value
                                                                     - 0.5F),
                                                                 0.F,
                                                                 1.F);
      colour[2] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * value, 0.F, 1.F);
    }

    float ColPixel::GetDensity() const
    {
      return density;
    }

    float ColPixel::GetStress() const
    {
      return stress;
    }

    unsigned int ColPixel::GetI() const
    {
      return i.i;
    }
    unsigned int ColPixel::GetJ() const
    {
      return i.j;
    }
    bool ColPixel::IsRT() const
    {
      return i.isRt;
    }
    bool ColPixel::IsGlyph() const
    {
      return i.isGlyph;
    }
    bool ColPixel::IsStreakline() const
    {
      return i.isStreakline;
    }

  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::ColPixel>::RegisterMpiDataType()
  {
    int col_pixel_count = 15;
    int col_pixel_blocklengths[15] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    MPI_Datatype col_pixel_types[15] = { MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_FLOAT,
                                         MPI_INT,
                                         MPI_INT,
                                         MPI_UB };

    MPI_Aint col_pixel_disps[15];

    col_pixel_disps[0] = 0;

    for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(float)
            * col_pixel_blocklengths[i - 1]);
      }
      else if (col_pixel_types[i - 1] == MPI_INT)
      {
        col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
      }
    }
    MPI_Datatype type;
    MPI_Type_struct(col_pixel_count,
                    col_pixel_blocklengths,
                    col_pixel_disps,
                    col_pixel_types,
                    &type);
    MPI_Type_commit(&type);
    return type;
  }
}
