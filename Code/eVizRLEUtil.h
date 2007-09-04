/*----------------------------------------------------------------------------

    (C)Copyright 2004 The University of Manchester, United Kingdom,
    all rights reserved.

    This software is produced by the Supercomputing, Visualization &
    e-Science Group, Manchester Computing, the Victoria University of
    Manchester as part of the e-Viz project.

    This software has been tested with care but is not guaranteed for
    any particular purpose. Neither the copyright holder, nor the
    University of Manchester offer any warranties or representations,
    nor do they accept any liabilities with respect to this software.

    This software must not be used for commercial gain without the
    written permission of the authors.

    This software must not be redistributed without the written
    permission of the authors.

    Permission is granted to modify this software, provided any
    modifications are made freely available to the original authors.

    Supercomputing, Visualization & e-Science Group
    Manchester Computing
    University of Manchester
    Manchester M13 9PL

    WWW:    http://www.sve.man.ac.uk
    email:  sve@man.ac.uk
    Tel:    +44 161 275 6095
    Fax:    +44 161 275 6800

    Initial version by: M Riding, 02.04.2004

---------------------------------------------------------------------------*/

#ifndef __EVIZRLEUTIL_H__
#define __EVIZRLEUTIL_H__

#ifdef __cplusplus
extern "C" {
#endif

// choose some arbitray number - not 0 or 255
#define ESCAPE_CHAR 42 


/** Compress an image in memory using the RLE codec. The function does not
 *  allocate any memory.
 *
 *  @param pixelData      pointer to the image to compress
 *  @param width          width of the image in pixels
 *  @param height         height of the image in pixels
 *  @param bpp            number of BYTES per pixel
 *  @param compressedSize the size of the compressed image in bytes
 *  @param compressedData return address of the compressed data (up to you to 
 *                        allocate a large enough buffer)
 *
 *  @return          an indication of success
 */
int eViz_RLE_writeMemory(const unsigned char *pixelData, 
			 const unsigned int width, const unsigned int height, 
			 const unsigned int bpp, 
			 int *compressedSize, unsigned char *compressedData);

/** Compress an image in memory using the RLED codec. The function does not
 *  allocate any memory.
 *
 *  @param pixelData      pointer to the image to compress
 *  @param previous       pointer to the uncompressed previous frame
 *  @param width          width of the image in pixels
 *  @param height         height of the image in pixels
 *  @param bpp            number of BYTES per pixel
 *  @param compressedSize the size of the compressed image in bytes
 *  @param compressedData return address of the compressed data (up to you to 
 *                        allocate a large enough buffer)
 *
 *  @return          an indication of success
 */
int eViz_RLED_writeMemory(const unsigned char *pixelData, 
			  const unsigned char *previous, 
			  const unsigned int width, 
			  const unsigned int height, 
			  const unsigned int bpp, 
			  int *compressedSize, unsigned char *compressedData);

/** Decompress a RLE image in memory. The function does not allocate any 
 *  memory.
 *
 *  @param compressedData pointer to the image to decompress
 *  @param compressedSize size of the image in memory in bytes
 *  @param width          width of the image in pixels
 *  @param height         height of the image in pixels
 *  @param bpp            number of BYTES per pixel
 *  @param pixelData      return address of the uncompressed data (up to you to
 *                        allocate a large enough buffer)
 *
 *  @return               an indication of success
 */
int eViz_RLE_readMemory(const unsigned char *compressedData, 
			const unsigned int compressedSize, 
			unsigned int *width, unsigned int *height, 
			unsigned int *bpp, unsigned char *pixelData);

/** Decompress a RLED image in memory. The function does not allocate any 
 *  memory.
 *
 *  @param compressedData pointer to the image to decompress
 *  @param previous       pointer to the uncompressed previous frame
 *  @param compressedSize size of the image in memory in bytes
 *  @param width          width of the image in pixels
 *  @param height         height of the image in pixels
 *  @param bpp            number of BYTES per pixel
 *  @param pixelData      return address of the uncompressed data (up to you to
 *                        allocate a large enough buffer)
 *
 *  @return               an indication of success
 */
int eViz_RLED_readMemory(const unsigned char *compressedData,
			 const unsigned char *previous,
			 const unsigned int compressedSize,
			 unsigned int *width, unsigned int *height, 
			 unsigned int *bpp, unsigned char *pixelData);

#ifdef __cplusplus
}
#endif

#endif
