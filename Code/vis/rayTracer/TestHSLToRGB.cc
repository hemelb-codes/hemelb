#include <iostream>

#include "HSLToRGBConverter.h"

int main( int argc, const char* argv[] )
{
  while (1)
  {
    float lHue;
    std::cout<< "Hue: " << std::endl;
    std::cin >> lHue;

    float lSaturation;
    std::cout << "Saturation: " << std::endl;
    std::cin >> lSaturation;

    float lLightness;
    std::cout << "Lightness: " << std::endl;
    std::cin >> lLightness;
 
    unsigned char lRGB[3];
    hemelb::vis::raytracer::HSLToRGBConverter::ConvertHSLToRGB
      (lHue, lSaturation, lLightness, lRGB);

    std::cout << "R: " << (int)lRGB[0] 
	      << " G: " << (int)lRGB[1] 
	      << " B: " << (int)lRGB[2] 
	      << std::endl;
 }

}
