#include <vector>
#include <math.h>

#include "util/utilityFunctions.h"
#include "vis/Control.h"
#include "vis/RayTracer.h"
#include "vis/GlyphDrawer.h"

#include "io/XdrFileWriter.h"

namespace hemelb
{
  namespace vis
  {
    Control::Control(lb::StressTypes iStressType, geometry::LatticeData* iLatDat)
    {
      mVisSettings.mStressType = iStressType;

      this->vis = new Vis;

      //sites_x etc are globals declared in net.h
      vis->half_dim[0] = 0.5F * float (iLatDat->GetXSiteCount());
      vis->half_dim[1] = 0.5F * float (iLatDat->GetYSiteCount());
      vis->half_dim[2] = 0.5F * float (iLatDat->GetZSiteCount());

      vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1], vis->half_dim[2]));

      mVisSettings.mouse_x = -1;
      mVisSettings.mouse_y = -1;
    }

    void Control::initLayers(topology::NetworkTopology * iNetworkTopology,
                             geometry::LatticeData* iLatDat)
    {
      myRayTracer = new RayTracer(iNetworkTopology, iLatDat, &mDomainStats, &mScreen, &mViewpoint,
                                  &mVisSettings);
      myGlypher = new GlyphDrawer(iLatDat, &mScreen, &mDomainStats, &mViewpoint, &mVisSettings);

#ifndef NO_STREAKLINES
      myStreaker = new StreaklineDrawer(iNetworkTopology, iLatDat, &mScreen, &mViewpoint,
                                        &mVisSettings);
#endif
      // Note that rtInit does stuff to this->ctr_x (because this has
      // to be global)
      mVisSettings.ctr_x -= vis->half_dim[0];
      mVisSettings.ctr_y -= vis->half_dim[1];
      mVisSettings.ctr_z -= vis->half_dim[2];
    }

    void Control::SetProjection(const int &iPixels_x,
                                const int &iPixels_y,
                                const float &iLocal_ctr_x,
                                const float &iLocal_ctr_y,
                                const float &iLocal_ctr_z,
                                const float &iLongitude,
                                const float &iLatitude,
                                const float &iZoom)
    {
      float rad = 5.F * vis->system_size;
      float dist = 0.5 * rad;

      float centre[3] = { iLocal_ctr_x, iLocal_ctr_y, iLocal_ctr_z };

      mViewpoint.SetViewpointPosition(iLongitude * DEG_TO_RAD, iLatitude * DEG_TO_RAD, centre, rad,
                                      dist);

      mScreen.Set( (0.5 * vis->system_size) / iZoom, (0.5 * vis->system_size) / iZoom, iPixels_x,
                  iPixels_y, rad, &mViewpoint);
    }

    void Control::RegisterSite(int i, float density, float velocity, float stress)
    {
      myRayTracer->UpdateClusterVoxel(i, density, velocity, stress);
    }

    void Control::SetSomeParams(const float iBrightness,
                                const float iDensityThresholdMin,
                                const float iDensityThresholdMinMaxInv,
                                const float iVelocityThresholdMaxInv,
                                const float iStressThresholdMaxInv)
    {
      mVisSettings.brightness = iBrightness;
      mDomainStats.density_threshold_min = iDensityThresholdMin;

      mDomainStats.density_threshold_minmax_inv = iDensityThresholdMinMaxInv;
      mDomainStats.velocity_threshold_max_inv = iVelocityThresholdMaxInv;
      mDomainStats.stress_threshold_max_inv = iStressThresholdMaxInv;
    }

    void Control::updateImageSize(int pixels_x, int pixels_y)
    {
      mScreen.Resize(pixels_x, pixels_y);
    }

    /**
     * Method that uses a binary tree communication method to composite the image
     * as created from multiple processors.
     *
     * @param recv_buffer_id
     * @param iNetTopology
     */
    void Control::compositeImage(int recv_buffer_id, const topology::NetworkTopology * iNetTopology)
    {
      mScreen.CompositeImage(recv_buffer_id, &mVisSettings, iNetTopology);
    }

    void Control::render(int recv_buffer_id,
                         geometry::LatticeData* iLatDat,
                         const topology::NetworkTopology* iNetTopology)
    {
      mScreen.Reset();

      myRayTracer->Render();

      if (mVisSettings.mode == VisSettings::ISOSURFACESANDGLYPHS)
      {
        myGlypher->Render();
      }
#ifndef NO_STREAKLINES
      if (mVisSettings.mStressType == lb::ShearStress || mVisSettings.mode
          == VisSettings::WALLANDSTREAKLINES)
      {
        myStreaker->render(iLatDat);
      }
#endif
      compositeImage(recv_buffer_id, iNetTopology);
    }

    void Control::writeImage(int recv_buffer_id, std::string image_file_name)
    {
      io::XdrFileWriter writer = io::XdrFileWriter(image_file_name);

      writer << (int) mVisSettings.mode;

      writer << mDomainStats.physical_pressure_threshold_min
          << mDomainStats.physical_pressure_threshold_max
          << mDomainStats.physical_velocity_threshold_max
          << mDomainStats.physical_stress_threshold_max;

      mScreen.WritePixelCount(recv_buffer_id, &writer);
      mScreen.WritePixels(recv_buffer_id, &mDomainStats, &mVisSettings, &writer);
    }

    void Control::setMouseParams(double iPhysicalPressure, double iPhysicalStress)
    {
      mVisSettings.mouse_pressure = iPhysicalPressure;
      mVisSettings.mouse_stress = iPhysicalStress;
    }

    bool Control::MouseIsOverPixel(float* density, float* stress)
    {
      if (mVisSettings.mouse_x < 0 || mVisSettings.mouse_y < 0)
      {
        return false;
      }

      return mScreen.MouseIsOverPixel(RECV_BUFFER_A, mVisSettings.mouse_x, mVisSettings.mouse_y,
                                      density, stress);
    }

    void Control::streaklines(int time_step, int period, geometry::LatticeData* iLatDat)
    {
      myStreaker ->StreakLines(time_step, period, iLatDat);
    }

    void Control::restart()
    {
      myStreaker->Restart();
    }

    Control::~Control()
    {
#ifndef NO_STREAKLINES
      delete myStreaker;
#endif

      delete vis;
      delete myGlypher;
      delete myRayTracer;
    }

  } // namespace vis
} // namespace hemelb
