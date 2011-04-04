#include "vis/Screen.h"

namespace hemelb
{
  namespace vis
  {
    /**
     * Add a pixel to the screen.
     *
     * @param newPixel The new pixel to be added
     * @param iStressType The stress type of the visualisation
     * @param mode Controls what aspects of the visualisation to display.
     */
    void Screen::AddPixel(const ColPixel* newPixel, lb::StressTypes iStressType, int mode)
    {
      // Get the id of the pixel if we've already added one at the same location.
      int pixelId = col_pixel_id[newPixel->i.i * PixelsY + newPixel->i.j];

      // If we have one at this location, merge in the pixel.
      if (pixelId != -1)
      {
        localPixels[pixelId].MergeIn(newPixel, iStressType, mode);
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
        col_pixel_id[newPixel->i.i * PixelsY + newPixel->i.j] = col_pixels;

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
                            lb::StressTypes iStressType,
                            int mode)
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
        RenderLineHelper<true> (x, y, dy, dy - dx, x2, iStressType, mode);
      }
      else
      {
        RenderLineHelper<false> (x, y, dx, dx - dy, y2, iStressType, mode);
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
                                  lb::StressTypes stressType,
                                  int mode)
    {
      int d = incNE;

      while ( (xLimited && x <= limit) || (!xLimited && y <= limit))
      {
        if (x >= 0 && x < PixelsX && y >= 0 && y < PixelsY)
        {
          ColPixel col_pixel;
          col_pixel.i = PixelId(x, y);
          col_pixel.i.isGlyph = true;

          AddPixel(&col_pixel, stressType, mode);
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
  }
}
