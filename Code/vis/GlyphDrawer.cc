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
      int n = -1;

      for (unsigned int i = 0; i < mLatDat->GetXSiteCount(); i += mLatDat->GetBlockSize())
      {
        for (unsigned int j = 0; j < mLatDat->GetYSiteCount(); j += mLatDat->GetBlockSize())
        {
          for (unsigned int k = 0; k < mLatDat->GetZSiteCount(); k += mLatDat->GetBlockSize())
          {
            n++;
            geometry::LatticeData::BlockData * map_block_p = mLatDat->GetBlock(n);

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            unsigned int site_i = (mLatDat->GetBlockSize() >> 1);
            unsigned int site_j = (mLatDat->GetBlockSize() >> 1);
            unsigned int site_k = (mLatDat->GetBlockSize() >> 1);

            unsigned int m = ( ( (site_i << mLatDat->GetLog2BlockSize()) + site_j)
                << mLatDat->GetLog2BlockSize()) + site_k;

            if (map_block_p->site_data[m] & BIG_NUMBER3)
            {
              continue;
            }

            Glyph *lGlyph = new Glyph();

            lGlyph->x = float (i + site_i) - 0.5F * float (mLatDat->GetXSiteCount());
            lGlyph->y = float (j + site_j) - 0.5F * float (mLatDat->GetYSiteCount());
            lGlyph->z = float (k + site_k) - 0.5F * float (mLatDat->GetZSiteCount());

            lGlyph->f = mLatDat->GetFOld(map_block_p->site_data[m] * D3Q15::NUMVECTORS);

            mGlyphs.push_back(lGlyph);

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

    // Function to perform the rendering.
    void GlyphDrawer::render()
    {
      float screen_max[4];
      screen_max[0] = mScreen->MaxXValue;
      screen_max[1] = mScreen->MaxXValue;
      screen_max[2] = mScreen->MaxYValue;
      screen_max[3] = mScreen->MaxYValue;

      float scale[4];
      scale[0] = scale[1] = mScreen->ScaleX;
      scale[2] = scale[3] = mScreen->ScaleY;

      double density;
      double vx, vy, vz;
      double temp;

      float p1[3], p2[3], p3[3], p4[3];

      for (unsigned int n = 0; n < mGlyphs.size(); n++)
      {
        D3Q15::CalculateDensityAndVelocity(mGlyphs[n]->f, density, vx, vy, vz);

        temp = glyph_length * mLatDat->GetBlockSize() * mDomainStats->velocity_threshold_max_inv
            / density;

        vx *= temp;
        vy *= temp;
        vz *= temp;

        p1[0] = mGlyphs[n]->x;
        p1[1] = mGlyphs[n]->y;
        p1[2] = mGlyphs[n]->z;

        p2[0] = mGlyphs[n]->x + vx;
        p2[1] = mGlyphs[n]->y + vy;
        p2[2] = mGlyphs[n]->z + vz;

        mViewpoint->Project(p1, p3);
        mViewpoint->Project(p2, p4);

        p3[0] = scale[0] * (p3[0] + screen_max[0]);
        p3[1] = scale[1] * (p3[1] + screen_max[1]);
        p4[0] = scale[2] * (p4[0] + screen_max[2]);
        p4[1] = scale[3] * (p4[1] + screen_max[3]);

        mScreen->RenderLine(p3, p4, mVisSettings->mStressType, mVisSettings->mode);
      }
    }

  }
}
