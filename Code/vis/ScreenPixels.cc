#include <string.h>
#include <stddef.h>

#include "debug/Debugger.h"
#include "vis/ScreenPixels.h"

namespace hemelb
{
  namespace vis
  {
    ScreenPixels::ScreenPixels()
    {
      PixelsX = 0;
      PixelsY = 0;

      Reset();
    }

    ScreenPixels::~ScreenPixels()
    {
    }

    void ScreenPixels::Reset()
    {
      pixelCount = 0;

      for (unsigned int ii = 0; ii < COLOURED_PIXELS_MAX; ++ii)
      {
        pixelId[ii] = -1;
      }
    }

    void ScreenPixels::FoldIn(const ScreenPixels* inScreen, const VisSettings* visSettings)
    {
      if (pixelCount == 0)
      {
        memcpy(pixels, inScreen->pixels, sizeof(ColPixel) * inScreen->pixelCount);
        pixelCount = inScreen->pixelCount;
        for (unsigned int ii = 0; ii < COLOURED_PIXELS_MAX; ++ii)
        {
          pixelId[ii] = inScreen->pixelId[ii];
        }
      }

      AddPixels(inScreen->pixels, inScreen->pixelCount, visSettings);
    }

    void ScreenPixels::AddPixel(const ColPixel* newPixel, const VisSettings* visSettings)
    {
      // Get the id of the pixel if we've already added one at the same location.
      int id = pixelId[newPixel->GetI() * PixelsY + newPixel->GetJ()];

      // If we have one at this location, merge in the pixel.
      if (id != -1)
      {
        pixels[id].MergeIn(newPixel, visSettings);
      }
      // Otherwise, if we have exceeded the maximum number of pixels, do nothing.

      else if (pixelCount >= COLOURED_PIXELS_MAX)
      {
        return;
      }
      // Otherwise, add the pixel to the list.

      else
      {
        // Put the pixel number into the store of ids.
        pixelId[newPixel->GetI() * PixelsY + newPixel->GetJ()] = pixelCount;

        // Add the pixel to the end of the list and move the end marker.
        pixels[pixelCount] = *newPixel;
        ++pixelCount;
      }
    }

    void ScreenPixels::AddPixels(const ColPixel* newPixel,
                                 unsigned int pixelCount,
                                 const VisSettings* visSettings)
    {
      for (unsigned int ii = 0; ii < pixelCount; ++ii)
      {
        AddPixel(&newPixel[ii], visSettings);
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
    void ScreenPixels::RenderLine(const float endPoint1[3],
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
    void ScreenPixels::RenderLineHelper(int x,
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

    void ScreenPixels::WriteImage(io::Writer* writer,
                                  const DomainStats* domainStats,
                                  const VisSettings* visSettings) const
    {
      *writer << (int) visSettings->mode;

      *writer << domainStats->physical_pressure_threshold_min
          << domainStats->physical_pressure_threshold_max
          << domainStats->physical_velocity_threshold_max
          << domainStats->physical_stress_threshold_max;

      *writer << GetPixelsX();
      *writer << GetPixelsY();
      *writer << pixelCount;

      WritePixels(writer, domainStats, visSettings);
    }

    void ScreenPixels::WritePixels(io::Writer* writer,
                                   const DomainStats* domainStats,
                                   const VisSettings* visSettings) const
    {
      const int bits_per_char = sizeof(char) * 8;

      for (unsigned int i = 0; i < pixelCount; i++)
      {
        const ColPixel col_pixel = pixels[i];
        // Use a ray-tracer function to get the necessary pixel data.
        int index;
        unsigned char rgb_data[12];
        col_pixel.rawWritePixel(&index, rgb_data, domainStats, visSettings);

        *writer << index;

        int pix_data[3];
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

    void ScreenPixels::SetSize(int x, int y)
    {
      if (x * y <= (int) COLOURED_PIXELS_MAX)
      {
        PixelsX = x;
        PixelsY = y;
      }

      for (unsigned int i = 0; i < PixelsX * PixelsY; i++)
      {
        pixelId[i] = -1;
      }
    }

    int ScreenPixels::GetPixelsX() const
    {
      return PixelsX;
    }
    int ScreenPixels::GetPixelsY() const
    {
      return PixelsY;
    }
  }
}
