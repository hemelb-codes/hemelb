#include <vector>
#include <math.h>
#include <limits>

#include "util/utilityFunctions.h"
#include "vis/Control.h"
#include "vis/RayTracer.h"
#include "vis/GlyphDrawer.h"

#include "io/XdrFileWriter.h"

namespace hemelb
{
  namespace vis
  {
    Control::Control(lb::StressTypes iStressType,
                     const topology::NetworkTopology* netTop,
                     geometry::LatticeData* iLatDat)
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

      initLayers(netTop, iLatDat);
    }

    void Control::initLayers(const topology::NetworkTopology * iNetworkTopology,
                             geometry::LatticeData* iLatDat)
    {
      // We don't have all the minima / maxima on one core, so we have to gather them.
      // NOTE this only happens once, during initialisation, otherwise it would be
      // totally unforgivable.
      unsigned int block_min_x = std::numeric_limits<unsigned int>::max();
      unsigned int block_min_y = std::numeric_limits<unsigned int>::max();
      unsigned int block_min_z = std::numeric_limits<unsigned int>::max();
      unsigned int block_max_x = std::numeric_limits<unsigned int>::min();
      unsigned int block_max_y = std::numeric_limits<unsigned int>::min();
      unsigned int block_max_z = std::numeric_limits<unsigned int>::min();

      int n = -1;

      for (unsigned int i = 0; i < iLatDat->GetXBlockCount(); i++)
      {
        for (unsigned int j = 0; j < iLatDat->GetYBlockCount(); j++)
        {
          for (unsigned int k = 0; k < iLatDat->GetZBlockCount(); k++)
          {
            n++;

            geometry::LatticeData::BlockData * lBlock = iLatDat->GetBlock((unsigned int) n);
            if (lBlock->ProcessorRankForEachBlockSite == NULL)
            {
              continue;
            }

            block_min_x = util::NumericalFunctions::min(block_min_x, i);
            block_min_y = util::NumericalFunctions::min(block_min_y, j);
            block_min_z = util::NumericalFunctions::min(block_min_z, k);
            block_max_x = util::NumericalFunctions::max(block_max_x, i);
            block_max_y = util::NumericalFunctions::max(block_max_y, j);
            block_max_z = util::NumericalFunctions::max(block_max_z, k);
          }
        }
      }

      unsigned int mins[3], maxes[3];
      unsigned int localMins[3], localMaxes[3];

      localMins[0] = block_min_x;
      localMins[1] = block_min_y;
      localMins[2] = block_min_z;

      localMaxes[0] = block_max_x;
      localMaxes[1] = block_max_y;
      localMaxes[2] = block_max_z;

      MPI_Allreduce(localMins, mins, 3, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, maxes, 3, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

      mVisSettings.ctr_x = 0.5F * iLatDat->GetBlockSize() * (mins[0] + maxes[0]);
      mVisSettings.ctr_y = 0.5F * iLatDat->GetBlockSize() * (mins[1] + maxes[1]);
      mVisSettings.ctr_z = 0.5F * iLatDat->GetBlockSize() * (mins[2] + maxes[2]);

      myRayTracer = new RayTracer(iNetworkTopology,
                                  iLatDat,
                                  &mDomainStats,
                                  &mScreen,
                                  &mViewpoint,
                                  &mVisSettings);
      myGlypher = new GlyphDrawer(iLatDat, &mScreen, &mDomainStats, &mViewpoint, &mVisSettings);

#ifndef NO_STREAKLINES
      myStreaker = new StreaklineDrawer(iNetworkTopology,
                                        iLatDat,
                                        &mScreen,
                                        &mViewpoint,
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

      mViewpoint.SetViewpointPosition(iLongitude * DEG_TO_RAD,
                                      iLatitude * DEG_TO_RAD,
                                      centre,
                                      rad,
                                      dist);

      mScreen.Set( (0.5 * vis->system_size) / iZoom,
                   (0.5 * vis->system_size) / iZoom,
                  iPixels_x,
                  iPixels_y,
                  rad,
                  &mViewpoint);
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

    void Control::UpdateImageSize(int pixels_x, int pixels_y)
    {
      mScreen.Resize(pixels_x, pixels_y);
    }

    void Control::Render(geometry::LatticeData* iLatDat,
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
      mScreen.CompositeImage(&mVisSettings, iNetTopology);
    }

    void Control::WriteImage(std::string image_file_name)
    {
      io::XdrFileWriter writer = io::XdrFileWriter(image_file_name);

      writer << (int) mVisSettings.mode;

      writer << mDomainStats.physical_pressure_threshold_min
          << mDomainStats.physical_pressure_threshold_max
          << mDomainStats.physical_velocity_threshold_max
          << mDomainStats.physical_stress_threshold_max;

      mScreen.WritePixelCount(&writer);
      mScreen.WritePixels(&mDomainStats, &mVisSettings, &writer);
    }

    void Control::SetMouseParams(double iPhysicalPressure, double iPhysicalStress)
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

      return mScreen.MouseIsOverPixel(mVisSettings.mouse_x, mVisSettings.mouse_y, density, stress);
    }

    void Control::ProgressStreaklines(int time_step, int period, geometry::LatticeData* iLatDat)
    {
      myStreaker ->StreakLines(time_step, period, iLatDat);
    }

    void Control::Reset()
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
