#ifndef HEMELB_VIS_GLYPHDRAWER_H
#define HEMELB_VIS_GLYPHDRAWER_H

#include "net.h"
#include "lb/LocalLatticeData.h"
#include "lb/GlobalLatticeData.h"
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
        GlyphDrawer(Net* net,
                    hemelb::lb::GlobalLatticeData* iGlobalLatDat,
                    hemelb::lb::LocalLatticeData* iLocalLatDat);
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

        std::vector<Glyph*> mGlyphs;
        hemelb::lb::GlobalLatticeData* mGlobalLatDat;
    };

  }
}

#endif //HEMELB_VIS_GLYPHDRAWER_H
