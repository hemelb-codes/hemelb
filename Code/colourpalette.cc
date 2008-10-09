#include <math.h>
#include "colourpalette.h"

void ColourPalette (float value, float col[])
{
  col[0] = fminf(1.F, value);
  col[1] = 0.;
  col[2] = fmaxf(0.F, 1.F - value);
}


