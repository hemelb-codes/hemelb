#include "vis/GlyphDrawer.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {

    double GlyphDrawer::glyph_length = -1.F;

    // Constructor
    GlyphDrawer::GlyphDrawer(geometry::LatticeData* iLatDat,
                             Screen* iScreen,
                             DomainStats* iDomainStats,
                             Viewpoint* iViewpoint,
                             VisSettings* iVisSettings) :
      mLatDat(iLatDat), mScreen(iScreen), mDomainStats(iDomainStats), mViewpoint(iViewpoint),
          mVisSettings(iVisSettings)
    {
      unsigned int n = 0;

      // Iterate over the first site in each block.
      for (unsigned int i = 0; i < mLatDat->GetXSiteCount(); i += mLatDat->GetBlockSize())
      {
        for (unsigned int j = 0; j < mLatDat->GetYSiteCount(); j += mLatDat->GetBlockSize())
        {
          for (unsigned int k = 0; k < mLatDat->GetZSiteCount(); k += mLatDat->GetBlockSize())
          {
            // Get the block data for this block - if it has no site data, move on.
            geometry::LatticeData::BlockData * map_block_p = mLatDat->GetBlock(n);

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            // We put the glyph at the site at the centre of the block...
            const unsigned int site_i = (mLatDat->GetBlockSize() >> 1);
            const unsigned int site_j = (mLatDat->GetBlockSize() >> 1);
            const unsigned int site_k = (mLatDat->GetBlockSize() >> 1);

            const unsigned int siteIdOnBlock = ( ( (site_i << mLatDat->GetLog2BlockSize()) + site_j)
                << mLatDat->GetLog2BlockSize()) + site_k;

            // ... (only if there's fluid there).
            if (map_block_p->site_data[siteIdOnBlock] & BIG_NUMBER3)
            {
              continue;
            }

            // Create a glyph at the desired location
            Glyph *lGlyph = new Glyph();

            lGlyph->x = float (i + site_i) - 0.5F * float (mLatDat->GetXSiteCount());
            lGlyph->y = float (j + site_j) - 0.5F * float (mLatDat->GetYSiteCount());
            lGlyph->z = float (k + site_k) - 0.5F * float (mLatDat->GetZSiteCount());

            lGlyph->f = mLatDat->GetFOld(map_block_p->site_data[siteIdOnBlock] * D3Q15::NUMVECTORS);

            mGlyphs.push_back(lGlyph);

            ++n;
          } // for k
        } // for j
      } // for i

    }

    // Destructor
    GlyphDrawer::~GlyphDrawer()
    {
      for (unsigned int ii = 0; ii < mGlyphs.size(); ii++)
      {
        delete mGlyphs[ii];
      }
    }

    /**
     * Perform the rendering for each glyph.
     */
    void GlyphDrawer::render()
    {
      // For each glyph...
      for (unsigned int n = 0; n < mGlyphs.size(); n++)
      {
        // ... get the density and velocity at that point...
        double density;
        double vx, vy, vz;
        D3Q15::CalculateDensityAndVelocity(mGlyphs[n]->f, density, vx, vy, vz);

        // ... calculate the velocity vector multiplier...
        const double temp = glyph_length * mLatDat->GetBlockSize()
            * mDomainStats->velocity_threshold_max_inv / density;

        // ... calculate the two ends of the line we're going to draw...
        float p1[3], p2[3];
        p1[0] = mGlyphs[n]->x;
        p1[1] = mGlyphs[n]->y;
        p1[2] = mGlyphs[n]->z;

        p2[0] = mGlyphs[n]->x + vx * temp;
        p2[1] = mGlyphs[n]->y + vy * temp;
        p2[2] = mGlyphs[n]->z + vz * temp;

        // ... transform to the location on the screen, and render.
        float p3[3], p4[3];
        mViewpoint->Project(p1, p3);
        mViewpoint->Project(p2, p4);

        mScreen->Transform<float> (p3, p3);
        mScreen->Transform<float> (p4, p4);

        mScreen->RenderLine(p3, p4, mVisSettings->mStressType, mVisSettings->mode);
      }
    }

  }
}
