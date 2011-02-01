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
    // TODO ye gods.
    // make a global controller
    Control *controller;

    Control::Control(lb::StressTypes iStressType, lb::GlobalLatticeData &iGlobLatDat)
    {
      mStressType = iStressType;

      this->vis = new Vis;

      //sites_x etc are globals declared in net.h
      vis->half_dim[0] = 0.5F * float (iGlobLatDat.GetXSiteCount());
      vis->half_dim[1] = 0.5F * float (iGlobLatDat.GetYSiteCount());
      vis->half_dim[2] = 0.5F * float (iGlobLatDat.GetZSiteCount());

      vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1], vis->half_dim[2]));
      col_pixels_max = COLOURED_PIXELS_MAX;

      col_pixel_recv[0] = new ColPixel[col_pixels_max];
      col_pixel_recv[1] = new ColPixel[col_pixels_max];

      pixels_max = COLOURED_PIXELS_MAX;
      col_pixel_id = new int[pixels_max];

      for (int i = 0; i < COLOURED_PIXELS_MAX; i++)
      {
        col_pixel_id[i] = -1;
      }

    }

    void Control::initLayers(topology::NetworkTopology * iNetworkTopology,
                             lb::GlobalLatticeData &iGlobLatDat,
                             lb::LocalLatticeData &iLocalLatDat)
    {
      myRayTracer = new RayTracer(iNetworkTopology, &iLocalLatDat, &iGlobLatDat);
      myGlypher = new GlyphDrawer(&iGlobLatDat, &iLocalLatDat);

#ifndef NO_STREAKLINES
      myStreaker = new StreaklineDrawer(iNetworkTopology, iLocalLatDat, iGlobLatDat);
#endif
      // Note that rtInit does stuff to this->ctr_x (because this has
      // to be global)
      this->ctr_x -= vis->half_dim[0];
      this->ctr_y -= vis->half_dim[1];
      this->ctr_z -= vis->half_dim[2];
    }

    void Control::RotateAboutXThenY(const float &iSinThetaX,
                                    const float &iCosThetaX,
                                    const float &iSinThetaY,
                                    const float &iCosThetaY,
                                    const float &iXIn,
                                    const float &iYIn,
                                    const float &iZIn,
                                    float &oXOut,
                                    float &oYOut,
                                    float &oZOut)
    {
      // A rotation of iThetaX clockwise about the x-axis
      // Followed by a rotation of iThetaY anticlockwise about the y-axis.
      // In matrices:
      //       (cos(iThetaY)  0 sin(iThetaY)) (1 0            0              )
      // Out = (0             1 0           ) (0 cos(-iThetaX) -sin(-iThetaX)) In
      //       (-sin(iThetaY) 0 cos(iThetaY)) (0 sin(-iThetaX) cos(-iThetaX) )
      //
      //       (Xcos(iThetaY) + Zsin(iThetaY)cos(iThetaX) - Ysin(iThetaY)sin(iThetaX))
      // Out = (Ycos(iThetaX) + Zsin(iThetaX)                                        )
      //       (Zcos(iThetaX)cos(iThetaY) - Ysin(iThetaX)cos(iThetaY) - Xsin(iThetaY))

      const float lTemp = iZIn * iCosThetaX - iYIn * iSinThetaX;

      oXOut = lTemp * iSinThetaY + iXIn * iCosThetaY;
      oYOut = iZIn * iSinThetaX + iYIn * iCosThetaX;
      oZOut = lTemp * iCosThetaY - iXIn * iSinThetaY;
    }

    void Control::UndoRotateAboutXThenY(const float &iSinThetaX,
                                        const float &iCosThetaX,
                                        const float &iSinThetaY,
                                        const float &iCosThetaY,
                                        const float &iXIn,
                                        const float &iYIn,
                                        const float &iZIn,
                                        float &oXOut,
                                        float &oYOut,
                                        float &oZOut)
    {
      // Just do using the other Rotate function, swapping for x and y
      // (to swap order) and giving the inverse of the sine values
      // (to swap direction).
      RotateAboutXThenY(-iSinThetaY, iCosThetaY, -iSinThetaX, iCosThetaX, iYIn, iXIn, iZIn, oYOut,
                        oXOut, oZOut);
    }

    void Control::project(float p1[], float p2[])
    {
      float x1[3], x2[3];

      for (int l = 0; l < 3; l++)
      {
        x1[l] = p1[l] - mViewpoint.x[l];
      }

      UndoRotateAboutXThenY(mViewpoint.SinXRotation, mViewpoint.CosXRotation,
                            mViewpoint.SinYRotation, mViewpoint.CosYRotation, x1[0], x1[1], x1[2],
                            x2[0], x2[1], x2[2]);

      float temp = mViewpoint.dist / (p2[2] = -x2[2]);

      p2[0] = temp * x2[0];
      p2[1] = temp * x2[1];
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

      RotateAboutXThenY(mViewpoint.SinXRotation, mViewpoint.CosXRotation, mViewpoint.SinYRotation,
                        mViewpoint.CosYRotation, mScreen.MaxXValue, 0.0F, 0.0F,
                        mScreen.UnitVectorProjectionX[0], mScreen.UnitVectorProjectionX[1],
                        mScreen.UnitVectorProjectionX[2]);

      RotateAboutXThenY(mViewpoint.SinXRotation, mViewpoint.CosXRotation, mViewpoint.SinYRotation,
                        mViewpoint.CosYRotation, 0.0F, mScreen.MaxYValue, 0.0F,
                        mScreen.UnitVectorProjectionY[0], mScreen.UnitVectorProjectionY[1],
                        mScreen.UnitVectorProjectionY[2]);

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

    /**
     * Merge data from the first ColPixel argument into the second
     * ColPixel argument.
     */
    void Control::MergePixels(const ColPixel *fromPixel, ColPixel *toPixel)
    {
      // Merge raytracing data

      if (fromPixel->i.isRt && toPixel->i.isRt)
      {
        // Both are ray-tracing
        toPixel->vel_r += fromPixel->vel_r;
        toPixel->vel_g += fromPixel->vel_g;
        toPixel->vel_b += fromPixel->vel_b;

        if (mStressType != lb::ShearStress)
        {
          toPixel->stress_r += fromPixel->stress_r;
          toPixel->stress_g += fromPixel->stress_g;
          toPixel->stress_b += fromPixel->stress_b;
        }

        toPixel->dt += fromPixel->dt;

        if (fromPixel->t < toPixel->t)
        {
          toPixel->t = fromPixel->t;
          toPixel->density = fromPixel->density;
          toPixel->stress = fromPixel->stress;
        }

      }
      else if (fromPixel->i.isRt && !toPixel->i.isRt)
      {
        // Only the 'from' merge-pixel is ray-tracing
        toPixel->vel_r = fromPixel->vel_r;
        toPixel->vel_g = fromPixel->vel_g;
        toPixel->vel_b = fromPixel->vel_b;

        if (mStressType != lb::ShearStress)
        {
          toPixel->stress_r = fromPixel->stress_r;
          toPixel->stress_g = fromPixel->stress_g;
          toPixel->stress_b = fromPixel->stress_b;
        }

        toPixel->t = fromPixel->t;
        toPixel->dt = fromPixel->dt;
        toPixel->density = fromPixel->density;
        toPixel->stress = fromPixel->stress;

        toPixel->i.isRt = true;
      }
      // Done merging ray-tracing - (last combinations would be if from-pixel has no ray-tracing data)

      // Now merge glyph data
      if (mStressType != lb::ShearStress && (mode == 0 || mode == 1))
      {
        if (fromPixel->i.isGlyph)
        {
          toPixel->i.isGlyph = true;
        }
      }
      else
      {
#ifndef NO_STREAKLINES
        // merge streakline data
        if (fromPixel->i.isStreakline)
        {
          if (!toPixel->i.isStreakline)
          {
            toPixel->particle_z = fromPixel->particle_z;
            toPixel->particle_vel = fromPixel->particle_vel;
            toPixel->particle_inlet_id = fromPixel->particle_inlet_id;

            toPixel->i.isStreakline = true;
          }
          else
          {
            if (fromPixel->particle_z < toPixel->particle_z)
            {
              toPixel->particle_z = fromPixel->particle_z;
              toPixel->particle_vel = fromPixel->particle_vel;
              toPixel->particle_inlet_id = fromPixel->particle_inlet_id;
            }
          }
        }
#endif
      }
    }

    void Control::RegisterSite(int i, float density, float velocity, float stress)
    {
      myRayTracer->UpdateClusterVoxel(i, density, velocity, stress);
    }

    void Control::writePixel(ColPixel *col_pixel_p)
    {
      int *col_pixel_id_p, i, j;

      i = col_pixel_p->i.i;
      j = col_pixel_p->i.j;

      col_pixel_id_p = &col_pixel_id[i * mScreen.PixelsY + j];

      if (*col_pixel_id_p != -1)
      {
        MergePixels(col_pixel_p, &col_pixel_send[*col_pixel_id_p]);

      }
      else
      { // col_pixel_id_p == -1

        if (col_pixels >= COLOURED_PIXELS_MAX)
        {
          return;
        }

        *col_pixel_id_p = col_pixels;

        memcpy(&col_pixel_send[col_pixels], col_pixel_p, sizeof(ColPixel));
        ++col_pixels;
      }

    }

    void Control::renderLine(float p1[], float p2[])
    {
      // Store end points of the line and 'current' point (x and y).
      int x1, y1, x2, y2;

      int x = int (p1[0]);
      int y = int (p1[1]);

      if (int (p2[0]) < int (p1[0]))
      {
        x1 = int (p2[0]);
        y1 = int (p2[1]);
        x2 = x;
        y2 = y;
      }
      else
      {
        x1 = x;
        y1 = y;

        x2 = int (p2[0]);
        y2 = int (p2[1]);
      }

      // Initialise dy with the absolute difference in y between endpoints of the line, and
      // m with the sign (-1 / +1) of the gradient.
      int dy, m;
      if (y1 < y2)
      {
        dy = y2 - y1;
        m = 1;

      }
      else
      {
        dy = y1 - y2;
        m = -1;
      }

      // Set dx with the difference between x-values at the endpoints.
      int dx = x2 - x1;

      // Set up the iteration in general terms.
      //int incE, d, incNE, whileVariable, whileLimit, otherVariable, otherVariableIncrement;

      if (dx > dy)
      {
        int incE = dy;
        int d = dy - dx;
        int incNE = d;

        while (x <= x2)
        {
          if (! (x < 0 || x >= mScreen.PixelsX || y < 0 || y >= mScreen.PixelsY))
          {
            ColPixel col_pixel;
            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            writePixel(&col_pixel);
          }

          if (d < 0)
          {
            d += incE;
          }
          else
          {
            d += incNE;
            y += m;
          }
          ++x;

        } // end while

      }
      else if (y1 < y2)
      {
        int incE = dx;
        int d = dx - dy;
        int incNE = d;

        while (y <= y2)
        {
          if (! (x < 0 || x >= mScreen.PixelsX || y < 0 || y >= mScreen.PixelsY))
          {
            ColPixel col_pixel;
            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            writePixel(&col_pixel);
          }

          if (d < 0)
          {
            d += incE;
          }
          else
          {
            d += incNE;
            ++x;
          }
          ++y;

        } // while

      }
      else
      {
        int incE = dx;
        int d = dx - dy;
        int incNE = d;

        while (y >= y2)
        {
          if (! (x < 0 || x >= mScreen.PixelsX || y < 0 || y >= mScreen.PixelsY))
          {
            ColPixel col_pixel;
            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            writePixel(&col_pixel);
          }

          if (d < 0)
          {
            d += incE;
          }
          else
          {
            d += incNE;
            ++x;
          }
          --y;

        } // while
      }
    }

    void Control::SetSomeParams(const float iBrightness,
                                const float iDensityThresholdMin,
                                const float iDensityThresholdMinMaxInv,
                                const float iVelocityThresholdMaxInv,
                                const float iStressThresholdMaxInv)
    {
      brightness = iBrightness;
      density_threshold_min = iDensityThresholdMin;

      density_threshold_minmax_inv = iDensityThresholdMinMaxInv;
      velocity_threshold_max_inv = iVelocityThresholdMaxInv;
      stress_threshold_max_inv = iStressThresholdMaxInv;
    }

    void Control::updateImageSize(int pixels_x, int pixels_y)
    {
      if (pixels_x * pixels_y > mScreen.PixelsX * mScreen.PixelsY)
      {
        pixels_max = pixels_x * pixels_y;
        col_pixel_id = (int *) realloc(col_pixel_id, sizeof(int) * pixels_max);
      }
      for (int i = 0; i < pixels_x * pixels_y; i++)
      {
        col_pixel_id[i] = -1;
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
      int comm_inc, send_id, recv_id;
      int i, j;
      int m, n;

      ColPixel *col_pixel1, *col_pixel2;

      if (!iNetTopology->IsCurrentProcTheIOProc())
      {
        memcpy(col_pixel_recv[recv_buffer_id], col_pixel_send, col_pixels * sizeof(ColPixel));
      }

      comm_inc = 1;
      m = 1;

      while (m < iNetTopology->GetProcessorCount())
      {
        m <<= 1;
        int start_id = 1;

        for (recv_id = start_id; recv_id < iNetTopology->GetProcessorCount();)
        {
          send_id = recv_id + comm_inc;

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
            MPI_Send(&col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);

            if (col_pixels > 0)
            {
              MPI_Send(col_pixel_send, col_pixels, ColPixel::getMpiType(), recv_id, 20,
                       MPI_COMM_WORLD);
            }

          }
          else
          {
            MPI_Recv(&col_pixels_temp, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD, &status);

            if (col_pixels_temp > 0)
            {
              MPI_Recv(col_pixel_send, col_pixels_temp, ColPixel::getMpiType(), send_id, 20,
                       MPI_COMM_WORLD, &status);
            }

            for (n = 0; n < col_pixels_temp; n++)
            {
              col_pixel1 = &col_pixel_send[n];
              i = col_pixel1->i.i;
              j = col_pixel1->i.j;

              if (* (col_pixel_id_p = &col_pixel_id[i * mScreen.PixelsY + j]) == -1)
              {
                col_pixel2 = &col_pixel_recv[recv_buffer_id][*col_pixel_id_p = col_pixels];

                memcpy(col_pixel2, col_pixel1, sizeof(ColPixel));
                ++col_pixels;

              }
              else
              {
                col_pixel2 = &col_pixel_recv[recv_buffer_id][*col_pixel_id_p];

                MergePixels(col_pixel1, col_pixel2);
              }

            }
          }
          if (m < iNetTopology->GetProcessorCount() && iNetTopology->GetLocalRank() == recv_id)
          {
            memcpy(col_pixel_send, col_pixel_recv[recv_buffer_id], col_pixels * sizeof(ColPixel));
          }

          recv_id += comm_inc << 1;
        }
        comm_inc <<= 1;
      }

      if (iNetTopology->GetLocalRank() == 1)
      {
        MPI_Send(&col_pixels, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);

        if (col_pixels > 0)
        {
          MPI_Send(col_pixel_recv[recv_buffer_id], col_pixels, ColPixel::getMpiType(), 0, 20,
                   MPI_COMM_WORLD);
        }

      }
      else if (iNetTopology->GetLocalRank() == 0)
      {
        MPI_Recv(&col_pixels, 1, MPI_INT, 1, 20, MPI_COMM_WORLD, &status);

        if (col_pixels > 0)
        {
          MPI_Recv(col_pixel_recv[recv_buffer_id], col_pixels, ColPixel::getMpiType(), 1, 20,
                   MPI_COMM_WORLD, &status);
        }

      }
    }

    void Control::render(int recv_buffer_id,
                         lb::GlobalLatticeData &iGlobLatDat,
                         const topology::NetworkTopology * iNetTopology)
    {
      if (mScreen.PixelsX * mScreen.PixelsY > pixels_max)
      {
        pixels_max = util::max(2 * pixels_max, mScreen.PixelsX * mScreen.PixelsY);

        col_pixel_id = (int *) realloc(col_pixel_id, sizeof(int) * pixels_max);
      }
      col_pixels = 0;

      myRayTracer->Render(mStressType);

      if (mode == 1)
      {
        myGlypher->render();
      }
#ifndef NO_STREAKLINES
      if (shouldDrawStreaklines && (mStressType == lb::ShearStress || mode == 2))
      {
        myStreaker->render(iGlobLatDat);
      }
#endif
      compositeImage(recv_buffer_id, iNetTopology);

      col_pixels_recv[recv_buffer_id] = col_pixels;

      for (int m = 0; m < col_pixels_recv[recv_buffer_id]; m++)
      {
        col_pixel_id[col_pixel_send[m].i.i * mScreen.PixelsY + col_pixel_send[m].i.j] = -1;
      }
    }

    void Control::writeImage(int recv_buffer_id,
                             std::string image_file_name,
                             void(*ColourPalette)(float value, float col[]))
    {
      io::XdrFileWriter writer = io::XdrFileWriter(image_file_name);

      writer << mode;

      writer << physical_pressure_threshold_min << physical_pressure_threshold_max
          << physical_velocity_threshold_max << physical_stress_threshold_max;

      writer << mScreen.PixelsX << mScreen.PixelsY << col_pixels_recv[recv_buffer_id];

      for (int n = 0; n < col_pixels_recv[recv_buffer_id]; n++)
      {
        writer.writePixel(&col_pixel_recv[recv_buffer_id][n], ColourPalette, mStressType);
      }
    }

    void Control::setMouseParams(double iPhysicalPressure, double iPhysicalStress)
    {
      mouse_pressure = iPhysicalPressure;
      mouse_stress = iPhysicalStress;
    }

    void Control::streaklines(int time_step,
                              int period,
                              lb::GlobalLatticeData &iGlobLatDat,
                              lb::LocalLatticeData &iLocalLatDat)
    {
      myStreaker ->StreakLines(time_step, period, iGlobLatDat, iLocalLatDat);
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

      delete[] col_pixel_id;
    }

  } // namespace vis
} // namespace hemelb
