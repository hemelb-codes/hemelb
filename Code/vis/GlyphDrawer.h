#ifndef HEMELB_VIS_GLYPHDRAWER_H
#define HEMELB_VIS_GLYPHDRAWER_H

#include "lb/LocalLatticeData.h"
#include "lb/GlobalLatticeData.h"

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

        // TODO there might be another way of increasing access to this. It
        // currently needs to be public because it is one of the steering
        // parameters.
        static double glyph_length;

        // Constructor and destructor
        GlyphDrawer(hemelb::lb::GlobalLatticeData* iGlobalLatDat,
                    hemelb::lb::LocalLatticeData* iLocalLatDat,
                    Screen* iScreen,
                    DomainStats* iDomainStats,
                    Viewpoint* iViewpoint,
                    VisSettings* iVisSettings);
        ~GlyphDrawer();

        // Function to perform the rendering.
        void render();

      private:
        // A struct to represent a single glyph.
        struct Glyph
        {
            float x, y, z;
            double *f;
        };

        Screen* mScreen;
        DomainStats* mDomainStats;
        Viewpoint* mViewpoint;
        VisSettings* mVisSettings;

        std::vector<Glyph*> mGlyphs;
        hemelb::lb::GlobalLatticeData* mGlobalLatDat;
    };

  }
}

#endif //HEMELB_VIS_GLYPHDRAWER_H
