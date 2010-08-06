#ifndef __vis_colourPalette_h_
#define __vis_colourPalette_h_

namespace vis {
  
  // Type for conversion functions.
  typedef void (ColourPaletteFunction)(float, float *);

  class ColourPalette {
  public:
    // Function to populate a RGB colour (col) with values dependant
    //on the value of t.
    //static void PickColour (float t, float col[]);
    static ColourPaletteFunction PickColour;
  };
  
}
#endif //__vis_colourPalette.h_
