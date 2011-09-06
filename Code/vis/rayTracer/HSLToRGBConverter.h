#ifndef HEMELB_VIS_HSLTORBGCONVERTER_H
#define HEMELB_VIS_HSLTORBGCONVERTER_H

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      
      class HSLToRGBConverter
      {
      public:
	//iHue is between 0.0F and 360.0F
	//the other two between 0.0F and 1.0F
	static void ConvertHSLToRGB
	  ( float iHue, 
	    float iSaturation, 
	    float iLightness,
	    unsigned char oRGBColour[3]);
      };
    }
  }
}

#endif // HEMELB_VIS_HSLTORBGCONVERTER_H
