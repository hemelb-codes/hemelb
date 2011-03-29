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

      col_pixel_recv[0] = new ColPixel[MAXCOLOUREDPIXELS];
      col_pixel_recv[1] = new ColPixel[MAXCOLOUREDPIXELS];

      pixels_max = MAXCOLOUREDPIXELS;
      mScreen.col_pixel_id = new int[pixels_max];

      for (int i = 0; i < MAXCOLOUREDPIXELS; i++)
      {
        mScreen.col_pixel_id[i] = -1;
      }

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
      mScreen.MaxXValue = (0.5 * vis->system_size) / iZoom;
      mScreen.MaxYValue = (0.5 * vis->system_size) / iZoom;

      mScreen.PixelsX = iPixels_x;
      mScreen.PixelsY = iPixels_y;

      // Convert to radians
      float temp = iLongitude * 0.01745329F;

      mViewpoint.SinYRotation = sinf(temp);
      mViewpoint.CosYRotation = cosf(temp);

      // Convert to radians
      temp = iLatitude * 0.01745329F;

      mViewpoint.SinXRotation = sinf(temp);
      mViewpoint.CosXRotation = cosf(temp);

      float rad = 5.F * vis->system_size;
      float dist = 0.5 * rad;

      temp = rad * mViewpoint.CosXRotation;

      mViewpoint.x[0] = temp * mViewpoint.SinYRotation + iLocal_ctr_x;
      mViewpoint.x[1] = rad * mViewpoint.SinXRotation + iLocal_ctr_y;
      mViewpoint.x[2] = temp * mViewpoint.CosYRotation + iLocal_ctr_z;

      mViewpoint.dist = dist;

      mViewpoint.RotateToViewpoint(mScreen.MaxXValue, 0.0F, 0.0F,
                                   &mScreen.UnitVectorProjectionX[0],
                                   &mScreen.UnitVectorProjectionX[1],
                                   &mScreen.UnitVectorProjectionX[2]);

      mViewpoint.RotateToViewpoint(0.0F, mScreen.MaxYValue, 0.0F,
                                   &mScreen.UnitVectorProjectionY[0],
                                   &mScreen.UnitVectorProjectionY[1],
                                   &mScreen.UnitVectorProjectionY[2]);

      mScreen.ScaleX = (float) iPixels_x / (2.F * mScreen.MaxXValue);
      mScreen.ScaleY = (float) iPixels_y / (2.F * mScreen.MaxYValue);

      temp = dist / rad;

      mScreen.vtx[0] = (mViewpoint.x[0] + temp * (iLocal_ctr_x - mViewpoint.x[0]))
          - mScreen.UnitVectorProjectionX[0] - mScreen.UnitVectorProjectionY[0] - mViewpoint.x[0];
      mScreen.vtx[1] = (mViewpoint.x[1] + temp * (iLocal_ctr_y - mViewpoint.x[1]))
          - mScreen.UnitVectorProjectionX[1] - mScreen.UnitVectorProjectionY[1] - mViewpoint.x[1];
      mScreen.vtx[2] = (mViewpoint.x[2] + temp * (iLocal_ctr_z - mViewpoint.x[2]))
          - mScreen.UnitVectorProjectionX[2] - mScreen.UnitVectorProjectionY[2] - mViewpoint.x[2];

      mScreen.UnitVectorProjectionX[0] *= (2.F / (float) iPixels_x);
      mScreen.UnitVectorProjectionX[1] *= (2.F / (float) iPixels_x);
      mScreen.UnitVectorProjectionX[2] *= (2.F / (float) iPixels_x);

      mScreen.UnitVectorProjectionY[0] *= (2.F / (float) iPixels_y);
      mScreen.UnitVectorProjectionY[1] *= (2.F / (float) iPixels_y);
      mScreen.UnitVectorProjectionY[2] *= (2.F / (float) iPixels_y);
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
      if (pixels_x * pixels_y > mScreen.PixelsX * mScreen.PixelsY)
      {
        pixels_max = pixels_x * pixels_y;
        mScreen.col_pixel_id = (int *) realloc(mScreen.col_pixel_id, sizeof(int) * pixels_max);
      }
      for (int i = 0; i < pixels_x * pixels_y; i++)
      {
        mScreen.col_pixel_id[i] = -1;
      }
    }

    void Control::compositeImage(int recv_buffer_id, const topology::NetworkTopology * iNetTopology)
    {
      // here, the communications needed to composite the image are
      // handled through a binary tree pattern and parallel pairwise
      // blocking communications.

      MPI_Status status;

      int *col_pixel_id_p;
      int col_pixels_temp;
      int i, j;
      int n;

      ColPixel *col_pixel1, *col_pixel2;

      if (!iNetTopology->IsCurrentProcTheIOProc())
      {
        memcpy(col_pixel_recv[recv_buffer_id], mScreen.localPixels, mScreen.col_pixels
            * sizeof(ColPixel));
      }

      unsigned int comm_inc = 1;
      unsigned int m = 1;

      while (m < iNetTopology->GetProcessorCount())
      {
        m <<= 1;
        unsigned int start_id = 1;

        for (unsigned int recv_id = start_id; recv_id < iNetTopology->GetProcessorCount();)
        {
          unsigned int send_id = recv_id + comm_inc;

          if (iNetTopology->GetLocalRank() != recv_id && iNetTopology->GetLocalRank() != send_id)
          {
            recv_id += comm_inc << 1;
            continue;
          }

          if (send_id >= iNetTopology->GetProcessorCount() || recv_id == send_id)
          {
            recv_id += comm_inc << 1;
            continue;
          }

          if (iNetTopology->GetLocalRank() == send_id)
          {
            MPI_Send(&mScreen.col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);

            if (mScreen.col_pixels > 0)
            {
              MPI_Send(mScreen.localPixels, mScreen.col_pixels, ColPixel::getMpiType(), recv_id,
                       20, MPI_COMM_WORLD);
            }

          }
          else
          {
            MPI_Recv(&col_pixels_temp, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD, &status);

            if (col_pixels_temp > 0)
            {
              MPI_Recv(mScreen.localPixels, col_pixels_temp, ColPixel::getMpiType(), send_id, 20,
                       MPI_COMM_WORLD, &status);
            }

            for (n = 0; n < col_pixels_temp; n++)
            {
              col_pixel1 = &mScreen.localPixels[n];
              i = col_pixel1->i.i;
              j = col_pixel1->i.j;

              if (* (col_pixel_id_p = &mScreen.col_pixel_id[i * mScreen.PixelsY + j]) == -1)
              {
                col_pixel2 = &col_pixel_recv[recv_buffer_id][*col_pixel_id_p = mScreen.col_pixels];

                memcpy(col_pixel2, col_pixel1, sizeof(ColPixel));
                ++mScreen.col_pixels;

              }
              else
              {
                col_pixel2 = &col_pixel_recv[recv_buffer_id][*col_pixel_id_p];

                col_pixel2->MergeIn(col_pixel1, mVisSettings.mStressType, mVisSettings.mode);
              }

            }
          }
          if (m < iNetTopology->GetProcessorCount() && iNetTopology->GetLocalRank() == recv_id)
          {
            memcpy(mScreen.localPixels, col_pixel_recv[recv_buffer_id], mScreen.col_pixels
                * sizeof(ColPixel));
          }

          recv_id += comm_inc << 1;
        }
        comm_inc <<= 1;
      }

      if (iNetTopology->GetLocalRank() == 1)
      {
        MPI_Send(&mScreen.col_pixels, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);

        if (mScreen.col_pixels > 0)
        {
          MPI_Send(col_pixel_recv[recv_buffer_id], mScreen.col_pixels, ColPixel::getMpiType(), 0,
                   20, MPI_COMM_WORLD);
        }

      }
      else if (iNetTopology->GetLocalRank() == 0)
      {
        MPI_Recv(&mScreen.col_pixels, 1, MPI_INT, 1, 20, MPI_COMM_WORLD, &status);

        if (mScreen.col_pixels > 0)
        {
          MPI_Recv(col_pixel_recv[recv_buffer_id], mScreen.col_pixels, ColPixel::getMpiType(), 1,
                   20, MPI_COMM_WORLD, &status);
        }

      }
    }

    void Control::render(int recv_buffer_id,
                         geometry::LatticeData* iLatDat,
                         const topology::NetworkTopology* iNetTopology)
    {
      if (mScreen.PixelsX * mScreen.PixelsY > pixels_max)
      {
        pixels_max = util::NumericalFunctions::max(2 * pixels_max, mScreen.PixelsX
            * mScreen.PixelsY);

        mScreen.col_pixel_id = (int *) realloc(mScreen.col_pixel_id, sizeof(int) * pixels_max);
      }
      mScreen.col_pixels = 0;

      myRayTracer->Render(mVisSettings.mStressType);

      if (mVisSettings.mode == 1)
      {
        myGlypher->render();
      }
#ifndef NO_STREAKLINES
      if (mVisSettings.mStressType == lb::ShearStress || mVisSettings.mode == 2)
      {
        myStreaker->render(iLatDat);
      }
#endif
      compositeImage(recv_buffer_id, iNetTopology);

      col_pixels_recv[recv_buffer_id] = mScreen.col_pixels;

      for (int m = 0; m < col_pixels_recv[recv_buffer_id]; m++)
      {
        mScreen.col_pixel_id[mScreen.localPixels[m].i.i * mScreen.PixelsY
            + mScreen.localPixels[m].i.j] = -1;
      }
    }

    void Control::writeImage(int recv_buffer_id,
                             std::string image_file_name,
                             void(*ColourPalette)(float value, float col[]))
    {
      io::XdrFileWriter writer = io::XdrFileWriter(image_file_name);

      writer << mVisSettings.mode;

      writer << mDomainStats.physical_pressure_threshold_min
          << mDomainStats.physical_pressure_threshold_max
          << mDomainStats.physical_velocity_threshold_max
          << mDomainStats.physical_stress_threshold_max;

      writer << mScreen.PixelsX << mScreen.PixelsY << col_pixels_recv[recv_buffer_id];

      for (int n = 0; n < col_pixels_recv[recv_buffer_id]; n++)
      {
        writer.writePixel(&col_pixel_recv[recv_buffer_id][n], ColourPalette, &mDomainStats,
                          mVisSettings.mode, mVisSettings.mStressType);
      }
    }

    void Control::setMouseParams(double iPhysicalPressure, double iPhysicalStress)
    {
      mVisSettings.mouse_pressure = iPhysicalPressure;
      mVisSettings.mouse_stress = iPhysicalStress;
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

      delete[] col_pixel_recv[0];
      delete[] col_pixel_recv[1];

      delete[] mScreen.col_pixel_id;
    }

  } // namespace vis
} // namespace hemelb
