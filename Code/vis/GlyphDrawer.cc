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
        mLatDat(iLatDat), mScreen(iScreen), mDomainStats(iDomainStats), mViewpoint(iViewpoint), mVisSettings(iVisSettings)
    {
      for (geometry::BlockTraverser blockTrav(*mLatDat); blockTrav.CurrentLocationValid(); blockTrav.TraverseOne())
      {
        // Get the block data for this block - if it has no site data, move on.
        const geometry::Block& block = blockTrav.GetCurrentBlockData();

        if (block.IsEmpty())
        {
          continue;
        }

        // We put the glyph at the site at the centre of the block...
        const util::Vector3D<site_t> midBlockSite = util::Vector3D<site_t>(mLatDat->GetBlockSize()) / (site_t) 2;

        const site_t siteIdOnBlock = mLatDat->GetLocalSiteIdFromLocalSiteCoords(midBlockSite);

        // ... (only if there's fluid there).
        if (block.SiteIsSolid(siteIdOnBlock))
        {
          continue;
        }

        // Create a glyph at the desired location
        Glyph lGlyph;

        util::Vector3D<site_t> globalSiteCoords = mLatDat->GetGlobalCoords(blockTrav.GetCurrentLocation(),
                                                                           midBlockSite);

        lGlyph.x = float(globalSiteCoords.x) - 0.5F * float(mLatDat->GetSiteDimensions().x);
        lGlyph.y = float(globalSiteCoords.y) - 0.5F * float(mLatDat->GetSiteDimensions().y);
        lGlyph.z = float(globalSiteCoords.z) - 0.5F * float(mLatDat->GetSiteDimensions().z);

        lGlyph.siteId = block.GetLocalContiguousIndexForSite(siteIdOnBlock);

        mGlyphs.push_back(lGlyph);
      }
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
        RenderLineHelper<true>(x, y, dy, dy - dx, x2, visSettings, pixelSet);
      }
      else
      {
        RenderLineHelper<false>(x, y, dx, dx - dy, y2, visSettings, pixelSet);
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
    PixelSet<BasicPixel>* GlyphDrawer::Render(const lb::MacroscopicPropertyCache& propertyCache)
    {
      PixelSet<BasicPixel>* pixelSet = GetUnusedPixelSet();

      pixelSet->Clear();

      // For each glyph...
      for (site_t n = 0; n < (site_t) mGlyphs.size(); n++)
      {
        // ... get the velocity at that point...
        const util::Vector3D<distribn_t>& velocity = propertyCache.velocityCache.Get(mGlyphs[n].siteId);

        // ... calculate the velocity vector multiplier...
        const double temp = mVisSettings->glyphLength * ((distribn_t) mLatDat->GetBlockSize())
            * mDomainStats->velocity_threshold_max_inv;

        // ... calculate the two ends of the line we're going to draw...
        util::Vector3D<float> p1 = util::Vector3D<float>(mGlyphs[n].x, mGlyphs[n].y, mGlyphs[n].z);
        util::Vector3D<float> p2 = p1
            + util::Vector3D<float>(float(velocity.x * temp), float(velocity.y * temp), float(velocity.z * temp));

        // ... transform to the location on the screen, and render.
        XYCoordinates<float> p3 = mViewpoint->FlatProject(p1);
        XYCoordinates<float> p4 = mViewpoint->FlatProject(p2);

        p3 = mScreen->TransformScreenToPixelCoordinates<float>(p3);
        p4 = mScreen->TransformScreenToPixelCoordinates<float>(p4);

        RenderLine(p3, p4, mVisSettings, pixelSet);
      }

      return pixelSet;
    }
  }
}
