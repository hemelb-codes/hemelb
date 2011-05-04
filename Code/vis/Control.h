#ifndef HEMELB_VIS_CONTROL_H
#define HEMELB_VIS_CONTROL_H

#include "geometry/LatticeData.h"
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
        Control(lb::StressTypes iStressType,
                const topology::NetworkTopology* netTop,
                geometry::LatticeData* iLatDat);
        ~Control();

        void SetSomeParams(const float iBrightness,
                           const distribn_t iDensityThresholdMin,
                           const distribn_t iDensityThresholdMinMaxInv,
                           const distribn_t iVelocityThresholdMaxInv,
                           const distribn_t iStressThresholdMaxInv);

        void SetProjection(const int &pixels_x,
                           const int &pixels_y,
                           const float &ctr_x,
                           const float &ctr_y,
                           const float &ctr_z,
                           const float &longitude,
                           const float &latitude,
                           const float &zoom);

        bool MouseIsOverPixel(float* density, float* stress);

        void ProgressStreaklines(unsigned long time_step, unsigned long period, geometry::LatticeData* iLatDat);
        void Reset();

        void UpdateImageSize(int pixels_x, int pixels_y);
        void Render(geometry::LatticeData* iLatDat, const topology::NetworkTopology* iNetTopology);
        void WriteImage(std::string image_file_name);
        void SetMouseParams(double iPhysicalPressure, double iPhysicalStress);
        void RegisterSite(site_t i, distribn_t density, distribn_t velocity, distribn_t stress);

        Screen mScreen;
        Viewpoint mViewpoint;
        DomainStats mDomainStats;
        VisSettings mVisSettings;

      private:
        struct Vis
        {
            float half_dim[3];
            float system_size;
        };

        void initLayers(const topology::NetworkTopology * iNetworkTopology,
                        geometry::LatticeData* iLatDat);

        Vis* vis;
        RayTracer *myRayTracer;
        GlyphDrawer *myGlypher;
        StreaklineDrawer *myStreaker;

    };
  }
}

#endif // HEMELB_VIS_CONTROL_H
