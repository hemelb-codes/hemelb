#include <stdlib.h>

#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include "vis/Screen.h"

namespace hemelb
{
  namespace vis
  {

    // TODO This is probably going to have to be cleverly redesigned. We need to pass the images around over several iterations without
    // interference between the steering and writte-to-disk images.

    Screen::Screen()
    {
      PixelsMax = COLOURED_PIXELS_MAX;
      col_pixel_id = new int[PixelsMax];

      for (unsigned int i = 0; i < COLOURED_PIXELS_MAX; i++)
      {
        col_pixel_id[i] = -1;
      }

      compositingBuffer = new ColPixel[COLOURED_PIXELS_MAX];
    }

    Screen::~Screen()
    {
      delete[] compositingBuffer;
      delete[] col_pixel_id;
    }

    /**
     * Add a pixel to the screen.
     *
     * @param newPixel The new pixel to be added
     * @param iStressType The stress type of the visualisation
     * @param mode Controls what aspects of the visualisation to display.
     */
    void Screen::AddPixel(const ColPixel* newPixel, const VisSettings* visSettings)
    {
      // Get the id of the pixel if we've already added one at the same location.
      int pixelId = col_pixel_id[newPixel->GetI() * PixelsY + newPixel->GetJ()];

      // If we have one at this location, merge in the pixel.
      if (pixelId != -1)
      {
        localPixels[pixelId].MergeIn(newPixel, visSettings);
      }
      // Otherwise, if we have exceeded the maximum number of pixels, do nothing.
      else if (col_pixels >= COLOURED_PIXELS_MAX)
      {
        return;
      }
      // Otherwise, add the pixel to the list.
      else
      {
        // Put the pixel number into the store of ids.
        col_pixel_id[newPixel->GetI() * PixelsY + newPixel->GetJ()] = col_pixels;

        // Add the pixel to the end of the list and move the end marker.
        localPixels[col_pixels] = *newPixel;
        ++col_pixels;
      }
    }

    /**
     * Render a line between two points on the screen.
     *
     * @param endPoint1
     * @param endPoint2
     * @param iStressType
     * @param mode
     */
    void Screen::RenderLine(const float endPoint1[3],
                            const float endPoint2[3],
                            const VisSettings* visSettings)
    {
      // Store end points of the line and 'current' point (x and y).
      int x = int (endPoint1[0]);
      int y = int (endPoint1[1]);

      int x1, y1, x2, y2;
      if (endPoint2[0] < endPoint1[0])
      {
        x1 = int (endPoint2[0]);
        y1 = int (endPoint2[1]);
        x2 = x;
        y2 = y;
      }
      else
      {
        x1 = x;
        y1 = y;
        x2 = int (endPoint2[0]);
        y2 = int (endPoint2[1]);
      }

      // Set dx with the difference between x-values at the endpoints.
      int dx = x2 - x1;

      // Set up the iteration in general terms.

      if (y2 <= y1)
      {
        int temp = y2;
        y2 = y;
        y = temp;
      }

      int dy = y2 - y;

      if (dx > dy)
      {
        RenderLineHelper<true> (x, y, dy, dy - dx, x2, visSettings);
      }
      else
      {
        RenderLineHelper<false> (x, y, dx, dx - dy, y2, visSettings);
      }
    }

    /**
     * Helper function for rendering a line. The template parameter indicates whether the end of
     * the line will be limited by the x or y dimension.
     *
     * @param x
     * @param y
     * @param incE
     * @param incNE
     * @param limit
     * @param stressType
     * @param mode
     */
    template<bool xLimited>
    void Screen::RenderLineHelper(int x,
                                  int y,
                                  int incE,
                                  int incNE,
                                  int limit,
                                  const VisSettings* visSettings)
    {
      int d = incNE;

      while ( (xLimited && x <= limit) || (!xLimited && y <= limit))
      {
        if (x >= 0 && x < (int) PixelsX && y >= 0 && y < (int) PixelsY)
        {
          ColPixel col_pixel(x, y);
          AddPixel(&col_pixel, visSettings);
        }

        if (d < 0)
        {
          d += incE;
          if (xLimited)
          {
            ++x;
          }
          else
          {
            ++y;
          }
        }
        else
        {
          d += incNE;
          ++y;
          ++x;
        }
      } // end while
    }

    void Screen::Set(float maxX,
                     float maxY,
                     int pixelsX,
                     int pixelsY,
                     float rad,
                     const Viewpoint* viewpoint)
    {
      MaxXValue = maxX;
      MaxYValue = maxX;

      viewpoint->RotateToViewpoint(MaxXValue, 0.0F, 0.0F, UnitVectorProjectionX);
      viewpoint->RotateToViewpoint(0.0F, MaxYValue, 0.0F, UnitVectorProjectionY);

      PixelsX = pixelsX;
      PixelsY = pixelsX;

      ScaleX = (float) pixelsX / (2.F * MaxXValue);
      ScaleY = (float) pixelsY / (2.F * MaxYValue);

      float radVector[3];
      viewpoint->RotateToViewpoint(0.F, 0.F, -rad, radVector);

      for (int ii = 0; ii < 3; ++ii)
      {
        vtx[ii] = (0.5F * radVector[ii]) - UnitVectorProjectionX[ii] - UnitVectorProjectionY[ii];
      }

      UnitVectorProjectionX[0] *= (2.F / (float) pixelsX);
      UnitVectorProjectionX[1] *= (2.F / (float) pixelsX);
      UnitVectorProjectionX[2] *= (2.F / (float) pixelsX);

      UnitVectorProjectionY[0] *= (2.F / (float) pixelsY);
      UnitVectorProjectionY[1] *= (2.F / (float) pixelsY);
      UnitVectorProjectionY[2] *= (2.F / (float) pixelsY);
    }

    void Screen::Resize(unsigned int newPixelsX, unsigned int newPixelsY)
    {
      if (newPixelsX * newPixelsY > PixelsX * PixelsY)
      {
        PixelsMax = util::NumericalFunctions::max(2 * PixelsMax, newPixelsX * newPixelsY);
        col_pixel_id = (int *) realloc(col_pixel_id, sizeof(int) * PixelsMax);
      }

      PixelsX = newPixelsX;
      PixelsY = newPixelsY;

      for (unsigned int i = 0; i < PixelsX * PixelsY; i++)
      {
        col_pixel_id[i] = -1;
      }
    }

    void Screen::Reset()
    {
      col_pixels = 0;
    }

    void Screen::CompositeImage(const VisSettings* visSettings)
    {
      // Status object for MPI comms.
      MPI_Status status;

      // For all processors with pixels, copy these to the receive buffer.
      for (unsigned int ii = 0; ii < col_pixels; ++ii)
      {
        compositingBuffer[ii] = localPixels[ii];
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

      topology::NetworkTopology* netTop = topology::NetworkTopology::Instance();

      // Start with a difference in rank of 1, doubling every time.
      for (proc_t deltaRank = 1; deltaRank < netTop->GetProcessorCount(); deltaRank <<= 1)
      {
        // The receiving proc is all the ranks that are 1 modulo (deltaRank * 2)
        for (proc_t receivingProc = 1; receivingProc < (netTop->GetProcessorCount() - deltaRank); receivingProc
            += deltaRank << 1)
        {
          proc_t sendingProc = receivingProc + deltaRank;

          // If we're the sending proc, do the send.
          if (netTop->GetLocalRank() == sendingProc)
          {
            MPI_Send(&col_pixels, 1, MpiDataType(col_pixels), receivingProc, 20, MPI_COMM_WORLD);

            if (col_pixels > 0)
            {
              MPI_Send(localPixels,
                       col_pixels,
                       MpiDataType<ColPixel> (),
                       receivingProc,
                       20,
                       MPI_COMM_WORLD);
            }
          }

          // If we're the receiving proc, receive.
          else if (netTop->GetLocalRank() == receivingProc)
          {
            unsigned int col_pixels_temp;

            MPI_Recv(&col_pixels_temp,
                     1,
                     MpiDataType(col_pixels_temp),
                     sendingProc,
                     20,
                     MPI_COMM_WORLD,
                     &status);

            if (col_pixels_temp > 0)
            {
              MPI_Recv(localPixels,
                       col_pixels_temp,
                       MpiDataType<ColPixel> (),
                       sendingProc,
                       20,
                       MPI_COMM_WORLD,
                       &status);
            }

            // Now merge the received pixels in with the local store of pixels.
            for (unsigned int n = 0; n < col_pixels_temp; n++)
            {
              ColPixel* col_pixel1 = &localPixels[n];

              int id = col_pixel1->GetI() * GetPixelsY() + col_pixel1->GetJ();
              if (col_pixel_id[id] == -1)
              {
                col_pixel_id[id] = col_pixels;

                compositingBuffer[col_pixels] = *col_pixel1;
                ++col_pixels;
              }
              else
              {
                compositingBuffer[col_pixel_id[id]].MergeIn(col_pixel1, visSettings);
              }
            }

            // If this isn't the last iteration, copy the pixels from the received buffer
            // back to the screen.
            if ( (deltaRank << 1) < netTop->GetProcessorCount())
            {
              for (unsigned int ii = 0; ii < col_pixels; ++ii)
              {
                localPixels[ii] = compositingBuffer[ii];
              }
            }
          }
        }
      }

      // Send the final image from proc 1 to 0.
      if (netTop->GetLocalRank() == 1)
      {
        MPI_Send(&col_pixels, 1, MpiDataType(col_pixels), 0, 20, MPI_COMM_WORLD);

        if (col_pixels > 0)
        {
          MPI_Send(compositingBuffer, col_pixels, MpiDataType<ColPixel> (), 0, 20, MPI_COMM_WORLD);
        }

      }
      // Receive the final image on proc 0.
      else if (netTop->GetLocalRank() == 0)
      {
        MPI_Recv(&col_pixels, 1, MpiDataType(col_pixels), 1, 20, MPI_COMM_WORLD, &status);

        if (col_pixels > 0)
        {
          MPI_Recv(compositingBuffer,
                   col_pixels,
                   MpiDataType<ColPixel> (),
                   1,
                   20,
                   MPI_COMM_WORLD,
                   &status);
        }
      }

      pixelCountInBuffer = col_pixels;

      for (unsigned int m = 0; m < pixelCountInBuffer; m++)
      {
        col_pixel_id[localPixels[m].GetI() * GetPixelsY() + localPixels[m].GetJ()] = -1;
      }
    }

    bool Screen::MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress)
    {
      for (unsigned int i = 0; i < pixelCountInBuffer; i++)
      {
        if (compositingBuffer[i].IsRT() && int (compositingBuffer[i].GetI()) == mouseX
            && int (compositingBuffer[i].GetJ()) == mouseY)
        {
          *density = compositingBuffer[i].GetDensity();
          *stress = compositingBuffer[i].GetStress();

          return true;
        }
      }

      return false;
    }

    void Screen::WritePixelCount(io::Writer* writer)
    {
      writer->operator <<(GetPixelsX());
      writer->operator <<(GetPixelsY());
      writer->operator <<(pixelCountInBuffer);
    }

    void Screen::WritePixels(const DomainStats* domainStats,
                             const VisSettings* visSettings,
                             io::Writer* writer)
    {
      int index;
      int pix_data[3];
      unsigned char rgb_data[12];
      int bits_per_char = sizeof(char) * 8;

      for (unsigned int i = 0; i < pixelCountInBuffer; i++)
      {
        //        col_pixel_p = &compositingBuffer[i];
        ColPixel& col_pixel = compositingBuffer[i];
        // Use a ray-tracer function to get the necessary pixel data.
        col_pixel.rawWritePixel(&index, rgb_data, domainStats, visSettings);

        *writer << index;

        pix_data[0] = (rgb_data[0] << (3 * bits_per_char)) + (rgb_data[1] << (2 * bits_per_char))
            + (rgb_data[2] << bits_per_char) + rgb_data[3];

        pix_data[1] = (rgb_data[4] << (3 * bits_per_char)) + (rgb_data[5] << (2 * bits_per_char))
            + (rgb_data[6] << bits_per_char) + rgb_data[7];

        pix_data[2] = (rgb_data[8] << (3 * bits_per_char)) + (rgb_data[9] << (2 * bits_per_char))
            + (rgb_data[10] << bits_per_char) + rgb_data[11];

        for (int i = 0; i < 3; i++)
        {
          *writer << pix_data[i];
        }
        *writer << io::Writer::eol;
      }
    }

    const float* Screen::GetVtx() const
    {
      return vtx;
    }
    const float* Screen::GetUnitVectorProjectionX() const
    {
      return UnitVectorProjectionX;
    }
    const float* Screen::GetUnitVectorProjectionY() const
    {
      return UnitVectorProjectionY;
    }
    int Screen::GetPixelsX() const
    {
      return PixelsX;
    }
    int Screen::GetPixelsY() const
    {
      return PixelsY;
    }
  }
}
