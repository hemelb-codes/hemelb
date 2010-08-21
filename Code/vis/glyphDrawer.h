#ifndef HEME_VIS_GLYPHDRAWER_H
#define HEME_VIS_GLYPHDRAWER_H

#include "vis/Control.h"
// TODO this could probably be reduced to the net class and some visualisation class.
#include "vis/rt.h"

namespace heme
{
  namespace vis
  {
    // Class for drawing glyphs.
    class GlyphDrawer : public Layer {
    public:
    
      // TODO there might be another way of increasing access to this. It
      // currently needs to be public because it is one of the steering
      // parameters.
      static double glyph_length;
  
      // Constructor and destructor
      GlyphDrawer(Net* net);
      ~GlyphDrawer();
  
      // Function to perform the rendering.
      virtual void render();
  
    private:
      // A struct to represent a single glyph.
      struct Glyph {
	float x, y, z;
	double *f;
      };
  
      // A pointer for an array of glyphs.
      Glyph *glyph;
      int glyphs;
  
    };

  }
}

#endif //HEME_VIS_GLYPHDRAWER_H
