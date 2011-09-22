#include <iostream>

#include "HSLToRGBConverter.h"

void TestColour( float iHue, 
		 float iSaturation, 
		 float iLightness,
		 float iExpectedRed,
		 float iExpectedGreen,
		 float iExpectedBlue )
{
  std::cout << "Testing HSL " 
	    << iHue << ", " 
	    << iSaturation << ", "
	    << iLightness << ": ";

  unsigned char lRGB[3];
  
  hemelb::vis::raytracer::HSLToRGBConverter::Convert
    (iHue, iSaturation, iLightness, lRGB);

  if (iExpectedRed == lRGB[0] &&
      iExpectedGreen == lRGB[1] &&
      iExpectedBlue == lRGB[2])
  {
    std::cout << "passed" << std::endl;
  }
  else
  {
    std::cout << "failed - got "
	      << static_cast<int>(lRGB[0]) << ", "
	      << static_cast<int>(lRGB[1]) << ", "
	      << static_cast<int>(lRGB[2]) << ", "
	      << std::endl;
  }
}

int main( int argc, const char* argv[] )
{
  TestColour(0, 0, 0, 0, 0, 0);
  TestColour(1, 1, 1, 255, 255, 255);
  
  TestColour(250, 0.5, 0.5, 84, 63, 191);
}
