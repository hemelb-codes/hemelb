#ifndef HEMELB_VIS_CONTROL_H
#define HEMELB_VIS_CONTROL_H

#include "net.h"

#include "vis/ColPixel.h"
#include "vis/GlyphDrawer.h"
#include "vis/StreaklineDrawer.h"
#include "vis/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    /* the last three digits of the pixel identifier are used to
     indicate if the pixel is coloured via the ray tracing technique
     and/or a glyph and/or a particle/pathlet */

    // Class to control and use the effects of different visualisation
    // methods.
    class Control
    {
      public:
        struct Screen
        {
            float vtx[3];

            // Projection of unit vectors along screen axes into normal space.
            float UnitVectorProjectionX[3];
            float UnitVectorProjectionY[3];

            float MaxXValue, MaxYValue;
            float ScaleX, ScaleY;

            int PixelsX, PixelsY;
        };

        struct Viewpoint
        {
            float x[3];
            float SinYRotation, CosYRotation;
            float SinXRotation, CosXRotation;
            float dist;
        };

        Control(float iStressType);
        ~Control();
        void initLayers(Net* net);

        void project(float p1[], float p2[]);
        void writePixel(ColPixel *col_pixel);
        void MergePixels(ColPixel *col_pixel1, ColPixel *col_pixel2);
        void RotateAboutXThenY(const float &iSinThetaX,
                               const float &iCosThetaX,
                               const float &iSinThetaY,
                               const float &iCosSinThetaY,
                               const float &iXIn,
                               const float &iYIn,
                               const float &iZIn,
                               float &oXOut,
                               float &oYOut,
                               float &oZOut);
        void UndoRotateAboutXThenY(const float &iSinThetaX,
                                   const float &iCosThetaX,
                                   const float &iSinThetaY,
                                   const float &iCosSinThetaY,
                                   const float &iXIn,
                                   const float &iYIn,
                                   const float &iZIn,
                                   float &oXOut,
                                   float &oYOut,
                                   float &oZOut);

        void SetSomeParams(const float iBrightness,
                           const float iDensityThresholdMin,
                           const float iDensityThresholdMinMaxInv,
                           const float iVelocityThresholdMaxInv,
                           const float iStressThresholdMaxInv);

        void SetProjection(const int &pixels_x,
                           const int &pixels_y,
                           const float &ctr_x,
                           const float &ctr_y,
                           const float &ctr_z,
                           const float &longitude,
                           const float &latitude,
                           const float &zoom);

        void renderLine(float x1[], float x2[]);

        void streaklines(int time_step, int period, Net *net);
        void restart();

        void updateImageSize(int pixels_x, int pixels_y);
        void render(int recv_buffer_id, Net *net);
        void writeImage(int recv_buffer_id,
                        char *image_file_name,
                        void(*ColourPalette)(float value, float col[]));
        void setMouseParams(double iPhysicalPressure, double iPhysicalStress);
        void compositeImage(int recv_buffer_id, Net *net);

        void RegisterSite(int i, float density, float velocity, float stress);

        // better public member vars than globals!
        int mode;
        int image_freq;
        int shouldDrawStreaklines;

        float velocity_threshold_max_inv;
        float stress_threshold_max_inv;
        float density_threshold_min, density_threshold_minmax_inv;
        float ctr_x, ctr_y, ctr_z;
        float streaklines_per_pulsatile_period, streakline_length;
        double mouse_pressure, mouse_stress;
        float brightness;

        Screen mScreen;
        Viewpoint mViewpoint;

        int col_pixels_recv[2]; // number received?
        ColPixel* col_pixel_recv[2];

        float physical_pressure_threshold_min;
        float physical_pressure_threshold_max;
        float physical_velocity_threshold_max;
        float physical_stress_threshold_max;

        int mouse_x, mouse_y;

      private:

        struct Vis
        {
            float half_dim[3];
            float system_size;
        };

        int col_pixels, // number of ColPixels (?)
            col_pixels_max; // max permitted of above

        int *col_pixel_id; // array of pixel IDs

        int time_step, cycle;
        int period, inlets;
        int pixels_max;

        int perform_rendering;

        ColPixel col_pixel_send[COLOURED_PIXELS_MAX];

        Vis* vis;

        RayTracer *myRayTracer;
        GlyphDrawer *myGlypher;
        StreaklineDrawer *myStreaker;

        float mStressType;

    };

    extern Control *controller;
  }
}

#endif // HEMELB_VIS_CONTROL_H
