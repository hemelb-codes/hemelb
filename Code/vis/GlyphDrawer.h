#ifndef HEMELB_VIS_GLYPHDRAWER_H
#define HEMELB_VIS_GLYPHDRAWER_H

#include "geometry/LatticeData.h"

#include "vis/Screen.h"
#include "vis/DomainStats.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

#include <vector>

namespace hemelb
{
  namespace vis
  {
    // Class for drawing glyphs.
    class GlyphDrawer
    {
      public:
        // Constructor and destructor
        GlyphDrawer(geometry::LatticeData* iLatDat,
                    Screen* iScreen,
                    DomainStats* iDomainStats,
                    Viewpoint* iViewpoint,
                    VisSettings* iVisSettings);
        ~GlyphDrawer();

        // Function to perform the rendering.
        void Render();

      private:
        // A struct to represent a single glyph.
        struct Glyph
        {
            float x, y, z;
            double *f;
        };

        geometry::LatticeData* mLatDat;

        Screen* mScreen;
        DomainStats* mDomainStats;
        Viewpoint* mViewpoint;
        VisSettings* mVisSettings;

        std::vector<Glyph*> mGlyphs;
    };

  }
}

#endif //HEMELB_VIS_GLYPHDRAWER_H
