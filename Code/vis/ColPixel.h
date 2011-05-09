#ifndef HEMELB_VIS_COLPIXEL_H
#define HEMELB_VIS_COLPIXEL_H

#include "mpiInclude.h"
#include "vis/DomainStats.h"
#include "vis/VisSettings.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace vis
  {
    class ColPixel
    {
      public:
        ColPixel();

        ColPixel(int i, int j);

        ColPixel(int i, int j, float particleVelocity, float particleZ, int particleInletId);

        ColPixel(int i,
                 int j,
                 float t,
                 float dt,
                 float density,
                 float stress,
                 const float velocityColour[3],
                 const float stressColour[3]);

        /**
         * Merge data from the first ColPixel argument into the second
         * ColPixel argument.
         */
        void MergeIn(const ColPixel* fromPixel, const VisSettings* visSettings);

        void rawWritePixel(int* pixel_index,
                           unsigned char rgb_data[12],
                           const DomainStats* iDomainStats,
                           const VisSettings* visSettings);

        float GetDensity() const;
        float GetStress() const;

        unsigned int GetI() const;
        unsigned int GetJ() const;

        bool IsRT() const;
        bool IsGlyph() const;
        bool IsStreakline() const;

        static void PickColour(float value, float colour[3]);

      private:
        struct PixelId
        {
            // Bitfield, fits in 4 bytes
            bool isRt :1;
            bool isGlyph :1;
            bool isStreakline :1;
            unsigned int i :14;
            unsigned int j :14;

            PixelId();
            PixelId(unsigned int i, unsigned int j);
        };

        void MakePixelColour(int rawRed, int rawGreen, int rawBlue, unsigned char* dest);

        // Pixel identity
        struct PixelId i;

        // Ray tracer pixel data
        float t, dt;
        float vel_r, vel_g, vel_b;
        float stress_r, stress_g, stress_b;
        float density;
        float stress;

        // Streakline pixel data
        float particle_vel;
        float particle_z;
        int particle_inlet_id;

    };

  }
}

#endif // HEMELB_VIS_COLPIXEL_H
