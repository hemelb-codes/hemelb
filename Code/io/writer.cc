#include "writer.h"
#include "../rt.h"

// Function to write out our struct, ColPixel.
void Writer::writePixel (ColPixel *col_pixel_p,
			 void (*colourPalette)(float value, float col[])) {
  // TODO: make this deal with spaces/newlines for general writer
  
  unsigned int index;
  unsigned int pix_data[3];
  unsigned char rgb_data[12];
  int bits_per_char = sizeof(char) * 8;
  
  // Use a ray-tracer function to get the necessary pixel data.
  rawWritePixel(col_pixel_p, &index, rgb_data, colourPalette);

  *this << index;
  
  pix_data[0] = (rgb_data[0] << (3*bits_per_char)) +
    (rgb_data[1] << (2*bits_per_char)) +
    (rgb_data[2] << bits_per_char) +
    rgb_data[3];
  
  pix_data[1] = (rgb_data[4] << (3*bits_per_char)) +
    (rgb_data[5] << (2*bits_per_char)) + 
    (rgb_data[6] << bits_per_char) +
    rgb_data[7];
  
  pix_data[2] = (rgb_data[8] << (3*bits_per_char)) + 
    (rgb_data[9] << (2*bits_per_char)) + 
    (rgb_data[10] << bits_per_char) + 
    rgb_data[11];
  
  for (int i=0; i<3; i++) {
    *this << pix_data[i];
  }
  *this << eol;
}


Writer& Writer::operator<< (enum Writer::Separator & value){
  writeRecordSeparator();
  return *this;
}
