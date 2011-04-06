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
    struct PixelId
    {
        unsigned int isRt :1;
        unsigned int isGlyph :1;
        unsigned int isStreakline :1;
        unsigned int i :14;
        unsigned int j :14;

        PixelId();
        PixelId(int i, int j);
    };

    class ColPixel
    {
      public:
        ColPixel();

        ColPixel(float particleVelocity, float particleZ, int particleInletId);

        ColPixel(float t,
                 float dt,
                 float density,
                 float stress,
                 const float velocityColour[3],
                 const float stressColour[3]);

        struct PixelId i;

        /**
         * Merge data from the first ColPixel argument into the second
         * ColPixel argument.
         */
        void
        MergeIn(const ColPixel* fromPixel, const VisSettings* visSettings);

        void rawWritePixel(int* pixel_index,
                           unsigned char rgb_data[12],
                           const DomainStats* iDomainStats,
                           const VisSettings* visSettings);
        static const MPI_Datatype& getMpiType();

        static void PickColour(float value, float colour[3]);

        float GetDensity();
        float GetStress();

      private:
        static void registerMpiType();
        static MPI_Datatype mpiType;

        void MakePixelColour(int rawRed, int rawGreen, int rawBlue, unsigned char* dest);

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
