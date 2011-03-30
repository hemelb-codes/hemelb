#ifndef HEMELB_VIS_CONTROL_H
#define HEMELB_VIS_CONTROL_H

#include "lb/GlobalLatticeData.h"
#include "lb/LocalLatticeData.h"
#include "lb/LbmParameters.h"

#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

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
        Control(lb::StressTypes iStressType, lb::GlobalLatticeData* iGlobLatDat);
        ~Control();

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

        void streaklines(int time_step,
                         int period,
                         lb::GlobalLatticeData* iGlobLatDat,
                         lb::LocalLatticeData* iLocalLatDat);
        void restart();

        void updateImageSize(int pixels_x, int pixels_y);
        void render(int recv_buffer_id,
                    lb::GlobalLatticeData* iGlobLatDat,
                    const topology::NetworkTopology* iNetTopology);
        void writeImage(int recv_buffer_id,
                        std::string image_file_name,
                        void(*ColourPalette)(float value, float col[]));
        void setMouseParams(double iPhysicalPressure, double iPhysicalStress);
        void compositeImage(int recv_buffer_id, const topology::NetworkTopology * iNetTopology);

        void RegisterSite(int i, float density, float velocity, float stress);

        void initLayers(topology::NetworkTopology * iNetworkTopology,
                        lb::GlobalLatticeData* iGlobLatDat,
                        lb::LocalLatticeData* iLocalLatDat);

        Screen mScreen;
        Viewpoint mViewpoint;
        DomainStats mDomainStats;
        VisSettings mVisSettings;

        int col_pixels_recv[2]; // number received?
        ColPixel* col_pixel_recv[2];

       private:
        struct Vis
        {
            float half_dim[3];
            float system_size;
        };

        int pixels_max;

        Vis* vis;

        RayTracer *myRayTracer;
        GlyphDrawer *myGlypher;
        StreaklineDrawer *myStreaker;

        static const long MAXCOLOUREDPIXELS = 2048 * 2048;

    };
  }
}

#endif // HEMELB_VIS_CONTROL_H
