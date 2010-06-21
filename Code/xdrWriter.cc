#include "xdrWriter.h"

// This constructor only exists so it can be set to protected to prevent instances of the base class being made.
// The alternative of declaring a pure virtual constructor or destructor doesn't work. Factory methods might be a 
// potential improvement here.
XdrWriter::XdrWriter()
{
}

void XdrWriter::writeInt(int* intToWrite)
{
  xdr_int(&myXdr, intToWrite);
}

void XdrWriter::writeDouble(double* doubleToWrite)
{
  xdr_double (&myXdr, doubleToWrite);
}

void XdrWriter::writeShort(short* shortToWrite)
{
  xdr_short (&myXdr, shortToWrite);
}

void XdrWriter::writeFloat(float* floatToWrite)
{
  xdr_float (&myXdr, floatToWrite);
}

void XdrWriter::writeUnsignedInt(unsigned int* unsignedIntToWrite)
{
  xdr_u_int (&myXdr, unsignedIntToWrite);
}

void XdrWriter::writePixel (ColPixel *col_pixel_p, void (*ColourPalette) (float value, float col[]))
{
  unsigned int index;
  unsigned int pix_data[3];
  unsigned char rgb_data[12];
  int bits_per_char = sizeof(char) * 8;
  
  rawWritePixel (col_pixel_p, &index, rgb_data, ColourPalette);

  writeUnsignedInt(&index);

  pix_data[0] = (rgb_data[0]<<(3*bits_per_char)) + (rgb_data[1]<<(2*bits_per_char)) + (rgb_data[2]<<bits_per_char) + rgb_data[3];
  pix_data[1] = (rgb_data[4]<<(3*bits_per_char)) + (rgb_data[5]<<(2*bits_per_char)) + (rgb_data[6]<<bits_per_char) + rgb_data[7];
  pix_data[2] = (rgb_data[8]<<(3*bits_per_char)) + (rgb_data[9]<<(2*bits_per_char)) + (rgb_data[10]<<bits_per_char) + rgb_data[11];

  writeUnsignedInt(&pix_data[0]);
  writeUnsignedInt(&pix_data[1]);
  writeUnsignedInt(&pix_data[2]);
}