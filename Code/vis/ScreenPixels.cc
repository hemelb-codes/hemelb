#include <string.h>
#include <stddef.h>
#include <stdlib.h>

#include "debug/Debugger.h"
#include "util/utilityFunctions.h"
#include "vis/ScreenPixels.h"

namespace hemelb
{
  namespace vis
  {
    ScreenPixels::ScreenPixels()
    {
      PixelsX = 0;
      PixelsY = 0;

      pixelArrayLength = 100;
      pixels = new ColPixel[pixelArrayLength];

      Reset();
    }

    ScreenPixels::~ScreenPixels()
    {
      delete[] pixels; 
    }

    void ScreenPixels::Reset()
    {
      storedPixels = 0;

      for (unsigned int ii = 0; ii < COLOURED_PIXELS_MAX; ++ii)
      {
        pixelId[ii] = -1;
      }
    }

    void ScreenPixels::FoldIn(const ScreenPixels* inScreen, const VisSettings* visSettings)
    {
      if (storedPixels == 0)
      {
        GuaranteeArraySize(inScreen->storedPixels);

        memcpy(pixels, inScreen->pixels, sizeof(ColPixel) * inScreen->storedPixels);
        storedPixels = inScreen->storedPixels;
        for (unsigned int ii = 0; ii < COLOURED_PIXELS_MAX; ++ii)
        {
          pixelId[ii] = inScreen->pixelId[ii];
        }
      }

      AddPixels(inScreen->pixels, inScreen->storedPixels, visSettings);
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

      else if (storedPixels >= COLOURED_PIXELS_MAX)
      {
        return;
      }
      // Otherwise, add the pixel to the list.

      else
      {
        // Put the pixel number into the store of ids.
        pixelId[newPixel->GetI() * PixelsY + newPixel->GetJ()] = storedPixels;

        GuaranteeArraySize(storedPixels + 1);

        // Add the pixel to the end of the list and move the end marker.
        pixels[storedPixels] = *newPixel;
        ++storedPixels;
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
    void ScreenPixels::RenderLine(const Vector3D<float>& endPoint1,
                                  const Vector3D<float>& endPoint2,
                                  const VisSettings* visSettings)
    {
      // Store end points of the line and 'current' point (x and y).
      int x = int(endPoint1.x);
      int y = int(endPoint1.y);

      int x1, y1, x2, y2;
      if (endPoint2.x < endPoint1.x)
      {
        x1 = int(endPoint2.x);
        y1 = int(endPoint2.y);
        x2 = x;
        y2 = y;
      }
      else
      {
        x1 = x;
        y1 = y;
        x2 = int(endPoint2.x);
        y2 = int(endPoint2.y);
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
      *writer << storedPixels;

      WritePixels(writer, domainStats, visSettings);
    }

    void ScreenPixels::WritePixels(io::Writer* writer,
                                   const DomainStats* domainStats,
                                   const VisSettings* visSettings) const
    {
      const int bits_per_char = sizeof(char) * 8;

      for (unsigned int i = 0; i < storedPixels; i++)
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

    void ScreenPixels::GuaranteeArraySize(unsigned int desiredSize)
    {
      if (pixelArrayLength < desiredSize)
      {
        // Expand 1.2 times the required size, to keep reallocs to a minimum.
        unsigned int newSize = desiredSize + desiredSize / 5;
        newSize = util::NumericalFunctions::min(newSize, COLOURED_PIXELS_MAX);

        pixels = (ColPixel*) realloc(pixels, sizeof(ColPixel) * newSize);
        pixelArrayLength = newSize;
      }
    }

    unsigned int ScreenPixels::GetStoredPixelCount() const
    {
      return storedPixels;
    }

    unsigned int* ScreenPixels::GetStoredPixelCountPtr()
    {
      return &storedPixels;
    }

    /**
     * Gets a pointer to the array of pixels.
     *
     * Because some functions will attempt to write to the array, but will have already
     * written the array length, guarantee that the array will be big enough for the incoming
     * count of stored pixels.
     * @return
     */
    ColPixel* ScreenPixels::GetPixelArray()
    {
      GuaranteeArraySize(GetStoredPixelCount());

      return pixels;
    }

  }
}
