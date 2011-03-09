#include "vis/Screen.h"

namespace hemelb
{
  namespace vis
  {
    void Screen::AddPixel(const ColPixel* newPixel, lb::StressTypes iStressType, int mode)
    {
      int *col_pixel_id_p = &col_pixel_id[newPixel->i.i * PixelsY + newPixel->i.j];

      if (*col_pixel_id_p != -1)
      {
        localPixels[*col_pixel_id_p].MergeIn(newPixel, iStressType, mode);
      }
      else
      { // col_pixel_id_p == -1

        if (col_pixels >= COLOURED_PIXELS_MAX)
        {
          return;
        }

        *col_pixel_id_p = col_pixels;

        memcpy(&localPixels[col_pixels], newPixel, sizeof(ColPixel));
        ++col_pixels;
      }
    }

    void Screen::RenderLine(const float endPoint1[3],
                            const float endPoint2[3],
                            lb::StressTypes iStressType,
                            int mode)
    {
      // Store end points of the line and 'current' point (x and y).
      int x1, y1, x2, y2;

      int x = int (endPoint1[0]);
      int y = int (endPoint1[1]);

      if (int (endPoint2[0]) < int (endPoint1[0]))
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
          if (! (x < 0 || x >= PixelsX || y < 0 || y >= PixelsY))
          {
            ColPixel col_pixel;
            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            AddPixel(&col_pixel, iStressType, mode);
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
          if (! (x < 0 || x >= PixelsX || y < 0 || y >= PixelsY))
          {
            ColPixel col_pixel;
            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            AddPixel(&col_pixel, iStressType, mode);
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
          if (! (x < 0 || x >= PixelsX || y < 0 || y >= PixelsY))
          {
            ColPixel col_pixel;
            col_pixel.i = PixelId(x, y);
            col_pixel.i.isGlyph = true;

            AddPixel(&col_pixel, iStressType, mode);
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

  }
}
