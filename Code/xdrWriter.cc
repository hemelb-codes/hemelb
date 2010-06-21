#include "xdrWriter.h"

XdrWriter::XdrWriter(char* fileName, char* failureString, char* successString)
{

}

void XdrWriter::xdrWritePixel (ColPixel *col_pixel_p, XDR* myXdr, void (*ColourPalette) (float value, float col[]))
{
  unsigned int index;
  unsigned int pix_data[3];
  unsigned char rgb_data[12];
  int bits_per_char = sizeof(char) * 8;
  
  rawWritePixel (col_pixel_p, &index, rgb_data, ColourPalette);

  xdr_u_int (myXdr, &index);

  pix_data[0] = (rgb_data[0]<<(3*bits_per_char)) + (rgb_data[1]<<(2*bits_per_char)) + (rgb_data[2]<<bits_per_char) + rgb_data[3];
  pix_data[1] = (rgb_data[4]<<(3*bits_per_char)) + (rgb_data[5]<<(2*bits_per_char)) + (rgb_data[6]<<bits_per_char) + rgb_data[7];
  pix_data[2] = (rgb_data[8]<<(3*bits_per_char)) + (rgb_data[9]<<(2*bits_per_char)) + (rgb_data[10]<<bits_per_char) + rgb_data[11];

  xdr_u_int (myXdr, &pix_data[0]);
  xdr_u_int (myXdr, &pix_data[1]);
  xdr_u_int (myXdr, &pix_data[2]);
}

XdrWriter::~XdrWriter()
{
 // xdr_destroy (&myXdr);
  fclose (myFile);
}