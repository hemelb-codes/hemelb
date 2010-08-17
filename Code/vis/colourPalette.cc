#include <math.h>
#include "vis/ColourPalette.h"

using namespace vis;

// Function to get an RGB colour in col, based on the value of the
// parameter t.  The colour is piece-wise linear in t.
void ColourPalette::PickColour (float t, float *col){
  if (t > 1.F) {
    col[0] = 1.F;
    col[1] = 0.F;
    col[2] = 0.F;
    
  } else if (t > 0.0F && t <= 0.25F) {
    col[0] = 0.F;
    col[1] = 4.F * t;
    col[2] = 1.F;
    
  } else if (t > 0.25F && t <= 0.5F) {
    col[0] = 0.F;
    col[1] = 1.F;
    col[2] = 2.F - 4.F * t;
    
  } else if (t > 0.5F && t <= 0.75F) {
    col[0] = 4.F * (t - 0.5F);
    col[1] = 1.F;
    col[2] = 0.F;
    
  } else if (t > 0.75F && t <= 1.0F) {
    col[0] = 1.F;
    col[1] = 4.F - 4.F * t;
    col[2] = 0.F;
  } else {
    col[0] = 0.F;
    col[1] = 0.F;
    col[2] = 1.F;
  }
  
}
