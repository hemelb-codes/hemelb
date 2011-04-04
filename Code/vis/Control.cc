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
      float rad = 5.F * vis->system_size;
      float dist = 0.5 * rad;

      float centre[3] = { iLocal_ctr_x, iLocal_ctr_y, iLocal_ctr_z };

      mViewpoint.SetViewpointPosition(iLongitude * DEG_TO_RAD, iLatitude * DEG_TO_RAD, centre, rad,
                                      dist);

      mViewpoint.RotateToViewpoint(mScreen.MaxXValue, 0.0F, 0.0F, mScreen.UnitVectorProjectionX);

      mViewpoint.RotateToViewpoint(0.0F, mScreen.MaxYValue, 0.0F, mScreen.UnitVectorProjectionY);

      mScreen.MaxXValue = (0.5 * vis->system_size) / iZoom;
      mScreen.MaxYValue = (0.5 * vis->system_size) / iZoom;

      mScreen.PixelsX = iPixels_x;
      mScreen.PixelsY = iPixels_y;

      mScreen.ScaleX = (float) iPixels_x / (2.F * mScreen.MaxXValue);
      mScreen.ScaleY = (float) iPixels_y / (2.F * mScreen.MaxYValue);

      float temp = dist / rad;

      const float* viewpointCentre = mViewpoint.GetViewpointCentre();

      mScreen.vtx[0] = (temp * (iLocal_ctr_x - viewpointCentre[0]))
          - mScreen.UnitVectorProjectionX[0] - mScreen.UnitVectorProjectionY[0];
      mScreen.vtx[1] = (temp * (iLocal_ctr_y - viewpointCentre[1]))
          - mScreen.UnitVectorProjectionX[1] - mScreen.UnitVectorProjectionY[1];
      mScreen.vtx[2] = (temp * (iLocal_ctr_z - viewpointCentre[2]))
          - mScreen.UnitVectorProjectionX[2] - mScreen.UnitVectorProjectionY[2];

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

    /**
     * Method that uses a binary tree communication method to composite the image
     * as created from multiple processors.
     *
     * @param recv_buffer_id
     * @param iNetTopology
     */
    void Control::compositeImage(int recv_buffer_id, const topology::NetworkTopology * iNetTopology)
    {
      // Status object for MPI comms.
      MPI_Status status;

      // For all processors with pixels, copy these to the receive buffer.
      if (!iNetTopology->IsCurrentProcTheIOProc())
      {
        for (unsigned int ii = 0; ii < mScreen.col_pixels; ++ii)
        {
          col_pixel_recv[recv_buffer_id][ii] = mScreen.localPixels[ii];
        }
      }

      /*
       * We do several iterations.
       *
       * On the first, every even proc passes data to the odd proc below, where it is merged.
       * On the second, the difference is two, so proc 3 passes to 1, 7 to 5, 11 to 9 etc.
       * On the third the differenec is four, so proc 5 passes to 1, 13 to 9 etc.
       * .
       * .
       * .
       *
       * This continues until all data is passed back to processor one, which passes it to proc 0.
       */

      // Start with a difference in rank of 1, doubling every time.
      for (unsigned int deltaRank = 1; deltaRank < iNetTopology->GetProcessorCount(); deltaRank
          <<= 1)
      {
        // The receiving proc is all the ranks that are 1 modulo (deltaRank * 2)
        for (unsigned int receivingProc = 1; receivingProc < (iNetTopology->GetProcessorCount()
            - deltaRank); receivingProc += deltaRank << 1)
        {
          unsigned int sendingProc = receivingProc + deltaRank;

          // If we're the sending proc, do the send.
          if (iNetTopology->GetLocalRank() == sendingProc)
          {
            MPI_Send(&mScreen.col_pixels, 1, MPI_UNSIGNED, receivingProc, 20, MPI_COMM_WORLD);

            if (mScreen.col_pixels > 0)
            {
              MPI_Send(mScreen.localPixels, mScreen.col_pixels, ColPixel::getMpiType(),
                       receivingProc, 20, MPI_COMM_WORLD);
            }
          }

          // If we're the receiving proc, receive.
          else if (iNetTopology->GetLocalRank() == receivingProc)
          {
            unsigned int col_pixels_temp;

            MPI_Recv(&col_pixels_temp, 1, MPI_UNSIGNED, sendingProc, 20, MPI_COMM_WORLD, &status);

            if (col_pixels_temp > 0)
            {
              MPI_Recv(mScreen.localPixels, col_pixels_temp, ColPixel::getMpiType(), sendingProc,
                       20, MPI_COMM_WORLD, &status);
            }

            // Now merge the received pixels in with the local store of pixels.
            for (unsigned int n = 0; n < col_pixels_temp; n++)
            {
              ColPixel* col_pixel1 = &mScreen.localPixels[n];

              int id = col_pixel1->i.i * mScreen.PixelsY + col_pixel1->i.j;
              if (mScreen.col_pixel_id[id] == -1)
              {
                mScreen.col_pixel_id[id] = mScreen.col_pixels;

                col_pixel_recv[recv_buffer_id][mScreen.col_pixels] = *col_pixel1;
                ++mScreen.col_pixels;
              }
              else
              {
                col_pixel_recv[recv_buffer_id][mScreen.col_pixel_id[id]].MergeIn(
                                                                                 col_pixel1,
                                                                                 mVisSettings.mStressType,
                                                                                 mVisSettings.mode);
              }
            }

            // If this isn't the last iteration, copy the pixels from the received buffer
            // back to the screen.
            if ( (deltaRank << 1) < iNetTopology->GetProcessorCount())
            {
              for (unsigned int ii = 0; ii < mScreen.col_pixels; ++ii)
              {
                mScreen.localPixels[ii] = col_pixel_recv[recv_buffer_id][ii];
              }
            }
          }
        }
      }

      // Send the final image from proc 1 to 0.
      if (iNetTopology->GetLocalRank() == 1)
      {
        MPI_Send(&mScreen.col_pixels, 1, MPI_UNSIGNED, 0, 20, MPI_COMM_WORLD);

        if (mScreen.col_pixels > 0)
        {
          MPI_Send(col_pixel_recv[recv_buffer_id], mScreen.col_pixels, ColPixel::getMpiType(), 0,
                   20, MPI_COMM_WORLD);
        }

      }
      // Receive the final image on proc 0.
      else if (iNetTopology->GetLocalRank() == 0)
      {
        MPI_Recv(&mScreen.col_pixels, 1, MPI_UNSIGNED, 1, 20, MPI_COMM_WORLD, &status);

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

      myRayTracer->Render();

      if (mVisSettings.mode == 1)
      {
        myGlypher->Render();
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
