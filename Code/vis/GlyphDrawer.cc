#include "vis/GlyphDrawer.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {

    double GlyphDrawer::glyph_length = -1.F;

    // Constructor
    GlyphDrawer::GlyphDrawer(hemelb::lb::GlobalLatticeData* iGlobalLatDat,
                             hemelb::lb::LocalLatticeData* iLocalLatDat)
    {
      mGlobalLatDat = iGlobalLatDat;
      int n = -1;

      for (unsigned int i = 0; i < iGlobalLatDat->GetXSiteCount(); i
          += iGlobalLatDat->GetBlockSize())
      {
        for (unsigned int j = 0; j < iGlobalLatDat->GetYSiteCount(); j
            += iGlobalLatDat->GetBlockSize())
        {
          for (unsigned int k = 0; k < iGlobalLatDat->GetZSiteCount(); k
              += iGlobalLatDat->GetBlockSize())
          {
            n++;
            lb::BlockData * map_block_p = &iGlobalLatDat->Blocks[n];

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            unsigned int site_i = (iGlobalLatDat->GetBlockSize() >> 1);
            unsigned int site_j = (iGlobalLatDat->GetBlockSize() >> 1);
            unsigned int site_k = (iGlobalLatDat->GetBlockSize() >> 1);

            unsigned int m = ( ( (site_i << iGlobalLatDat->Log2BlockSize) + site_j)
                << iGlobalLatDat->Log2BlockSize) + site_k;

            if (map_block_p->site_data[m] & (1U << 31U))
            {
              continue;
            }

            Glyph *lGlyph = new Glyph();

            lGlyph->x = float (i + site_i) - 0.5F * float (iGlobalLatDat->GetXSiteCount());
            lGlyph->y = float (j + site_j) - 0.5F * float (iGlobalLatDat->GetYSiteCount());
            lGlyph->z = float (k + site_k) - 0.5F * float (iGlobalLatDat->GetZSiteCount());

            int c1Plusc2 = 15;

            lGlyph->f = &iLocalLatDat->FOld[map_block_p->site_data[m] * c1Plusc2];

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
      screen_max[0] = vis::controller->mScreen.MaxXValue;
      screen_max[1] = vis::controller->mScreen.MaxXValue;
      screen_max[2] = vis::controller->mScreen.MaxYValue;
      screen_max[3] = vis::controller->mScreen.MaxYValue;

      float scale[4];
      scale[0] = scale[1] = vis::controller->mScreen.ScaleX;
      scale[2] = scale[3] = vis::controller->mScreen.ScaleY;

      double density;
      double vx, vy, vz;
      double temp;

      float p1[3], p2[3], p3[3], p4[3];

      for (unsigned int n = 0; n < mGlyphs.size(); n++)
      {
        D3Q15::CalculateDensityAndVelocity(mGlyphs[n]->f, density, vx, vy, vz);

        temp = glyph_length * mGlobalLatDat->GetBlockSize()
            * vis::controller->velocity_threshold_max_inv / density;

        vx *= temp;
        vy *= temp;
        vz *= temp;

        p1[0] = mGlyphs[n]->x;
        p1[1] = mGlyphs[n]->y;
        p1[2] = mGlyphs[n]->z;

        p2[0] = mGlyphs[n]->x + vx;
        p2[1] = mGlyphs[n]->y + vy;
        p2[2] = mGlyphs[n]->z + vz;

        vis::controller->project(p1, p3);
        vis::controller->project(p2, p4);

        p3[0] = scale[0] * (p3[0] + screen_max[0]);
        p3[1] = scale[1] * (p3[1] + screen_max[1]);
        p4[0] = scale[2] * (p4[0] + screen_max[2]);
        p4[1] = scale[3] * (p4[1] + screen_max[3]);

        vis::controller->renderLine(p3, p4);
      }
    }

  }
}
