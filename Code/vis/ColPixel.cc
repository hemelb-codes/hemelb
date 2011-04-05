#include <math.h>

#include "vis/ColPixel.h"
#include "util/utilityFunctions.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {

    MPI_Datatype MPI_col_pixel_type;
    PixelId::PixelId(int i_, int j_) :
      isRt(false), isGlyph(false), isStreakline(false), i(i_), j(j_)
    {
    }

    PixelId::PixelId() :
      isRt(false), isGlyph(false), isStreakline(false), i(0), j(0)
    {
    }

    MPI_Datatype ColPixel::mpiType = MPI_DATATYPE_NULL;

    // create the derived datatype for the MPI communications
    void ColPixel::registerMpiType()
    {
      int col_pixel_count = 15;
      int col_pixel_blocklengths[15] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

      MPI_Datatype col_pixel_types[15] = { MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
                                           MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
                                           MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_INT, MPI_UB };

      MPI_Aint col_pixel_disps[15];

      col_pixel_disps[0] = 0;

      for (int i = 1; i < col_pixel_count; i++)
      {
        if (col_pixel_types[i - 1] == MPI_FLOAT)
        {
          col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(float) * col_pixel_blocklengths[i
              - 1]);
        }
        else if (col_pixel_types[i - 1] == MPI_INT)
        {
          col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int)
              * col_pixel_blocklengths[i - 1]);
        }
      }
      MPI_Type_struct(col_pixel_count, col_pixel_blocklengths, col_pixel_disps, col_pixel_types,
                      &mpiType);
      MPI_Type_commit(&mpiType);
    }

    const MPI_Datatype& ColPixel::getMpiType()
    {
      if (mpiType == MPI_DATATYPE_NULL)
      {
        registerMpiType();
      }
      return mpiType;
    }

    void ColPixel::makePixelColour(unsigned char& red,
                                   unsigned char& green,
                                   unsigned char& blue,
                                   int rawRed,
                                   int rawGreen,
                                   int rawBlue)
    {
      red = (unsigned char) util::NumericalFunctions::enforceBounds(rawRed, 0, 255);
      green = (unsigned char) util::NumericalFunctions::enforceBounds(rawGreen, 0, 255);
      blue = (unsigned char) util::NumericalFunctions::enforceBounds(rawBlue, 0, 255);
    }

    /**
     * Merge data from the ColPixel argument into this pixel.
     */
    void ColPixel::MergeIn(const ColPixel *fromPixel, lb::StressTypes iStressType, int mode)
    {
      // Merge raytracing data

      if (fromPixel->i.isRt && i.isRt)
      {
        // Both are ray-tracing
        vel_r += fromPixel->vel_r;
        vel_g += fromPixel->vel_g;
        vel_b += fromPixel->vel_b;

        if (iStressType != lb::ShearStress)
        {
          stress_r += fromPixel->stress_r;
          stress_g += fromPixel->stress_g;
          stress_b += fromPixel->stress_b;
        }

        dt += fromPixel->dt;

        if (fromPixel->t < t)
        {
          t = fromPixel->t;
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

        if (iStressType != lb::ShearStress)
        {
          stress_r = fromPixel->stress_r;
          stress_g = fromPixel->stress_g;
          stress_b = fromPixel->stress_b;
        }

        t = fromPixel->t;
        dt = fromPixel->dt;
        density = fromPixel->density;
        stress = fromPixel->stress;

        i.isRt = true;
      }
      // Done merging ray-tracing - (last combinations would be if from-pixel has no ray-tracing data)

      // Now merge glyph data
      if (iStressType != lb::ShearStress && (mode == 0 || mode == 1))
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
                                 int mode,
                                 unsigned char rgb_data[],
                                 const DomainStats* iDomainStats,
                                 lb::StressTypes iLbmStressType)
    {
      float density_col[3], stress_col[3], particle_col[3];

      int bits_per_char = sizeof(char) * 8;
      int pixel_i, pixel_j;

      unsigned char r1, g1, b1;
      unsigned char r2, g2, b2;
      unsigned char r3, g3, b3;
      unsigned char r4, g4, b4;

      // store pixel id
      pixel_i = i.i;
      pixel_j = i.j;

      *pixel_index = (pixel_i << (2 * bits_per_char)) + pixel_j;

      r1 = g1 = b1 = 255;
      r2 = g2 = b2 = 255;

      if (i.isRt)
      {
        // store velocity volume rendering colour
        makePixelColour(r1, g1, b1, int (vel_r / dt), int (vel_g / dt), int (vel_b / dt));

        if (iLbmStressType != lb::ShearStress)
        {
          // store von Mises stress volume rendering colour
          makePixelColour(r2, g2, b2, int (stress_r / dt), int (stress_g / dt), int (stress_b / dt));

        }
        else if (stress < ((float) BIG_NUMBER))
        {
          PickColour(stress, stress_col);

          // store wall shear stress colour
          makePixelColour(r2, g2, b2, int (255.0F * stress_col[0]), int (255.0F * stress_col[1]),
                          int (255.0F * stress_col[2]));

        }
        else
        {
          r2 = g2 = b2 = 0;
        }

      } // if (isRt)

      if (iLbmStressType != lb::ShearStress && mode == 0)
      {
        PickColour(density, density_col);
        PickColour(stress, stress_col);

        // store wall pressure colour
        makePixelColour(r3, g3, b3, int (255.0F * density_col[0]), int (255.0F * density_col[1]),
                        int (255.0F * density_col[2]));

        // store von Mises stress colour
        makePixelColour(r4, g4, b4, int (255.0F * stress_col[0]), int (255.0F * stress_col[1]),
                        int (255.0F * stress_col[2]));

      }
      else if (iLbmStressType != lb::ShearStress && mode == 1)
      {
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
          makePixelColour(r3, g3, b3, int (127.5F * density_col[0]), int (127.5F * density_col[1]),
                          int (127.5F * density_col[2]));

          // store von Mises stress (+glyph) colour
          makePixelColour(r4, g4, b4, int (127.5F * stress_col[0]), int (127.5F * stress_col[1]),
                          int (127.5F * stress_col[2]));
        }
        else
        {
          r3 = g3 = b3 = 0;
          r4 = g4 = b4 = 0;
        }

      }
      else
      {

        if (i.isStreakline)
        {
          float scaled_vel = particle_vel * iDomainStats->velocity_threshold_max_inv;

          PickColour(scaled_vel, particle_col);

          // store particle colour
          makePixelColour(r3, g3, b3, int (255.0F * particle_col[0]),
                          int (255.0F * particle_col[1]), int (255.0F * particle_col[2]));

          r4 = r3;
          g4 = g3;
          b4 = b3;

        }
        else
        {
          // store pressure colour
          r3 = g3 = b3 = (unsigned char) util::NumericalFunctions::enforceBounds(int (127.5F
              * density), 0, 127);

          // store shear stress or von Mises stress
          if (stress < ((float) BIG_NUMBER))
          {
            r4 = g4 = b4 = (unsigned char) util::NumericalFunctions::enforceBounds(int (127.5F
                * stress), 0, 127);

          }
          else
          {
            r4 = g4 = b4 = 0;
          }
        }

      }

      rgb_data[0] = r1;
      rgb_data[1] = g1;
      rgb_data[2] = b1;
      rgb_data[3] = r2;
      rgb_data[4] = g2;
      rgb_data[5] = b2;
      rgb_data[6] = r3;
      rgb_data[7] = g3;
      rgb_data[8] = b3;
      rgb_data[9] = r4;
      rgb_data[10] = g4;
      rgb_data[11] = b4;

    }

    void ColPixel::PickColour(float value, float colour[3])
    {
      colour[0] = util::NumericalFunctions::enforceBounds<float>(4.F * value - 2.F, 0.F, 1.F);
      colour[1] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * fabs(value - 0.5F),
                                                                 0.F, 1.F);
      colour[2] = util::NumericalFunctions::enforceBounds<float>(2.F - 4.F * value, 0.F, 1.F);
    }

  }
}
