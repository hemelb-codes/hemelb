#include "vis/GlyphDrawer.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {
    // Constructor
    GlyphDrawer::GlyphDrawer(geometry::LatticeData* iLatDat,
                             Screen* iScreen,
                             DomainStats* iDomainStats,
                             Viewpoint* iViewpoint,
                             VisSettings* iVisSettings) :
      mLatDat(iLatDat), mScreen(iScreen), mDomainStats(iDomainStats), mViewpoint(iViewpoint),
          mVisSettings(iVisSettings)
    {
      int n = -1;

      // Iterate over the first site in each block.
      for (site_t i = 0; i < mLatDat->GetXSiteCount(); i += mLatDat->GetBlockSize())
      {
        for (site_t j = 0; j < mLatDat->GetYSiteCount(); j += mLatDat->GetBlockSize())
        {
          for (site_t k = 0; k < mLatDat->GetZSiteCount(); k += mLatDat->GetBlockSize())
          {
            ++n;

            // Get the block data for this block - if it has no site data, move on.
            geometry::BlockData * map_block_p = mLatDat->GetBlock(n);

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            // We put the glyph at the site at the centre of the block...
            const site_t site_i = (mLatDat->GetBlockSize() >> 1);
            const site_t site_j = (mLatDat->GetBlockSize() >> 1);
            const site_t site_k = (mLatDat->GetBlockSize() >> 1);

            const site_t siteIdOnBlock = ( ( (site_i << mLatDat->GetLog2BlockSize()) + site_j)
                << mLatDat->GetLog2BlockSize()) + site_k;

            // ... (only if there's fluid there).
            if (map_block_p->site_data[siteIdOnBlock] & BIG_NUMBER3)
            {
              continue;
            }

            // Create a glyph at the desired location
            Glyph lGlyph;

            lGlyph.x = float(i + site_i) - 0.5F * float(mLatDat->GetXSiteCount());
            lGlyph.y = float(j + site_j) - 0.5F * float(mLatDat->GetYSiteCount());
            lGlyph.z = float(k + site_k) - 0.5F * float(mLatDat->GetZSiteCount());

            lGlyph.f = mLatDat->GetFOld(map_block_p->site_data[siteIdOnBlock] * D3Q15::NUMVECTORS);

            mGlyphs.push_back(lGlyph);
          } // for k
        } // for j
      } // for i

    }

    /**
     * Destructor
     */
    GlyphDrawer::~GlyphDrawer()
    {
    }

    /**
     * Render a line between two points on the screen.
     *
     * @param endPoint1
     * @param endPoint2
     * @param iStressType
     * @param mode
     */
    void GlyphDrawer::RenderLine(const XYCoordinates<float>& endPoint1,
                                 const XYCoordinates<float>& endPoint2,
                                 const VisSettings* visSettings,
                                 PixelSet<BasicPixel>* pixelSet)
    {
      // Store end points of the line and 'current' point (x and y).
      int x = int(endPoint1.x);
      int y = int(endPoint1.y);

      // Ensure increasing x from point 1 to point 2.
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

      // Ensure increasing y.
      if (y2 <= y1)
      {
        int temp = y2;
        y2 = y;
        y = temp;
      }

      // Set dx with the difference between x-values at the endpoints.
      int dx = x2 - x1;
      int dy = y2 - y;

      if (dx > dy)
      {
        RenderLineHelper<true> (x, y, dy, dy - dx, x2, visSettings, pixelSet);
      }
      else
      {
        RenderLineHelper<false> (x, y, dx, dx - dy, y2, visSettings, pixelSet);
      }
    }

    /**
     * Helper function for rendering a line. The template parameter indicates whether the end of
     * the line will be limited by the x or y dimension.
     *
     * @tparam xLimited True if the line is limited by an x value, false if limited by a y value.
     * @param x Starting, minimal value of x.
     * @param y Starting, minimal value of y.
     * @param smallerDelta is the smaller of the two deltas between initial- and end- x and y.
     * @param deltaDifference is (smaller delta) - (larger delta)
     * @param limit Maximal value in the limiting direction
     */
    template<bool xLimited>
    void GlyphDrawer::RenderLineHelper(int x,
                                       int y,
                                       const int smallerDelta,
                                       const int deltaDifference,
                                       const int limit,
                                       const VisSettings* visSettings,
                                       PixelSet<BasicPixel>* pixelSet)
    {
      // This variable tracks whether we are above or below the line we draw.
      int overOrUnderMeasure = deltaDifference;

      while ( (xLimited && x <= limit) || (!xLimited && y <= limit))
      {
        // If on screen, add a pixel.
        if (x >= 0 && x < (int) mScreen->GetPixelsX() && y >= 0 && y < (int) mScreen->GetPixelsY())
        {
          BasicPixel pixel(x, y);
          pixelSet->AddPixel(pixel);
        }

        // We are effectively adding the smallerDelta every time (and increment in the other direction),
        // and we periodically subtract a largerDelta (and increment in the limited direction).
        if (overOrUnderMeasure < 0)
        {
          overOrUnderMeasure += smallerDelta;
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
          overOrUnderMeasure += deltaDifference;
          ++y;
          ++x;
        }
      } // end while
    }

    /**
     * Perform the rendering for each glyph.
     */
    PixelSet<BasicPixel>* GlyphDrawer::Render()
    {
      PixelSet<BasicPixel>* pixelSet = GetUnusedPixelSet();

      pixelSet->Clear();

      // For each glyph...
      for (site_t n = 0; n < (site_t) mGlyphs.size(); n++)
      {
        // ... get the density and velocity at that point...
        distribn_t density;
        distribn_t vx, vy, vz;
        D3Q15::CalculateDensityAndVelocity(mGlyphs[n].f, density, vx, vy, vz);

        // ... calculate the velocity vector multiplier...
        const double temp = mVisSettings->glyphLength * ((distribn_t) mLatDat->GetBlockSize())
            * mDomainStats->velocity_threshold_max_inv / density;

        // ... calculate the two ends of the line we're going to draw...
        util::Vector3D<float> p1 = util::Vector3D<float>(mGlyphs[n].x, mGlyphs[n].y, mGlyphs[n].z);
        util::Vector3D<float> p2 = p1 + util::Vector3D<float>(vx * temp, vy * temp, vz * temp);

        // ... transform to the location on the screen, and render.
        XYCoordinates<float> p3 = mViewpoint->FlatProject(p1);
        XYCoordinates<float> p4 = mViewpoint->FlatProject(p2);

        p3 = mScreen->TransformScreenToPixelCoordinates<float> (p3);
        p4 = mScreen->TransformScreenToPixelCoordinates<float> (p4);

        RenderLine(p3, p4, mVisSettings, pixelSet);
      }

      return pixelSet;
    }
  }
}
