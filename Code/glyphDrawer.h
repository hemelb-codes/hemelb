#ifndef __glyphDrawer_h_
#define __glyphDrawer_h_

#include "visualisationControl.h"
// TODO this could probably be reduced to the net class and some visualisation class.
#include "rt.h"

class glyphDrawer : public visualisationLayer
{
  private:

    struct Glyph
    {
      float x, y, z;
      double *f;
    };

    Glyph *glyph;
    int glyphs;

  public:
    static double vis_glyph_length;

    glyphDrawer(Net* net);
    ~glyphDrawer();
    virtual void render();
};

#endif //__glyphDrawer_h_