#ifndef __glyphDrawer_h_
#define __glyphDrawer_h_

#include "visualisationControl.h"
// TODO this could probably be reduced to the net class and some visualisation class.
#include "rt.h"

// Class for drawing glyphs.
class glyphDrawer : public visualisationLayer
{
  private:
    // A struct to represent a single glyph.
    struct Glyph
    {
      float x, y, z;
      double *f;
    };

    // A pointer for an array of glyphs.
    Glyph *glyph;
    int glyphs;

  public:
    // TODO there might be another way of increasing access to this. It currently needs to be public because it
    // is one of the steering parameters.
    static double vis_glyph_length;

    // Constructor and destructor
    glyphDrawer(Net* net);
    ~glyphDrawer();
    
    // Function to perform the rendering.
    virtual void render();
};

#endif //__glyphDrawer_h_