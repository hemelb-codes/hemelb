
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

// #include "eVizPNGUtil.h"
#include "eVizRLEUtil.h"

#include "eVizTypes.h"

#include <stdio.h>
#include <stdlib.h>

/* ***************************************************************** */
/* Local function definitions                                        */
/* ***************************************************************** */

/* int peek(int i, unsigned char *array, int arraySize); */

/* ***************************************************************** */
/* Externally accessible functions                                   */
/* ***************************************************************** */

int eViz_RLE_writeMemory(const unsigned char *pixelData,
                        const unsigned int width, const unsigned int height,
                        const unsigned int bpp,
                        int *compressedSize, unsigned char *compressedData)
{
 int s=0, d=0; // counters for the IO data
 int sSize = width * height * bpp;

 int sameByteCount = 1;
 unsigned char t,t2;
 unsigned short shortWidth;
 unsigned short shortHeight;
 unsigned char charBpp;

 int overFlowFlag = 0;

 if (pixelData == NULL || compressedSize == NULL || compressedData == NULL)
 {
   fprintf(stderr, "eViz_RLE_writeMemory: one or more arguments was NULL\n");
   return EVIZ_ERROR;
 }

 // first of all stick in the image dimensions
 // short, short,  byte = width, height, bpp
 shortWidth = width;
 shortHeight = height;
 charBpp = bpp;

 compressedData[0] = ((unsigned char*)&shortWidth)[0];
 compressedData[1] = ((unsigned char*)&shortWidth)[1];
 compressedData[2] = ((unsigned char*)&shortHeight)[0];
 compressedData[3] = ((unsigned char*)&shortHeight)[1];
 compressedData[4] = charBpp;

 d = 5;

 for (s=0; s<sSize; s++){
   int topt;

   // read bytes in, counting the number of identical consecutive bytes
   t = pixelData[s];

   // peak at the next bytes until they're different or we reach the end of
   // the data stream
   //t2 = peek(s+sameByteCount, pixelData, sSize);
   topt = s+sameByteCount;
   if (topt>=sSize)
   {
     t2 = -1;
     overFlowFlag = 1;
   }
   else
     t2 = pixelData[topt];

   while (t2 == t && !overFlowFlag && (sameByteCount < 254)){
     sameByteCount++;
     //t2 = peek(s+sameByteCount, pixelData, sSize);
     topt = s+sameByteCount;
     if (topt>=sSize)
     {
       t2 = -1;
 overFlowFlag = 1;
     }
     else
       t2 = pixelData[topt];
   }

   // encode data to stream

   if (sameByteCount == 1){
     // check to see if the byte is the escape char
     if (t == ESCAPE_CHAR){
       compressedData[d] = ESCAPE_CHAR;
       compressedData[d+1] = 1;
       compressedData[d+2] = t;
       d+=3;
     }
     else {
       compressedData[d] = t;
       d++;
     }
   }
   else if (sameByteCount == 2 && t != ESCAPE_CHAR){
     compressedData[d] = t;
     compressedData[d+1] = t;
     d+=2;
     s+=1; // (sameByteCount - 1)
   }
   else {
     compressedData[d] = ESCAPE_CHAR;
     compressedData[d+1] = sameByteCount;
     compressedData[d+2] = t;
     d+=3;
     s+=sameByteCount - 1;
   }

   sameByteCount = 1;

 }

 *compressedSize = d;

 return EVIZ_SUCCESS;
}

int eViz_RLED_writeMemory(const unsigned char *pixelData,
                         const unsigned char *previous,
                         const unsigned int width,
                         const unsigned int height,
                         const unsigned int bpp,
                         int *compressedSize, unsigned char *compressedData)
{
 int s=0, d=0; // counters for the IO data
 int sSize = width * height * bpp;

 int sameByteCount = 1;
 short diff, diff2;
 //unsigned char currentFrameTemp, previousFrameTemp;
 int topt;
 unsigned short shortWidth;
 unsigned short shortHeight;
 unsigned char charBpp;

 if (pixelData == NULL || previous == NULL ||
     compressedSize == NULL || compressedData == NULL)
 {
   fprintf(stderr, "eViz_RLED_writeMemory: one or more arguments was NULL\n");
   return EVIZ_ERROR;
 }

 // first of all stick in the image dimensions
 // short, short,  byte = width, height, bpp
 shortWidth = width;
 shortHeight = height;
 charBpp = bpp;

 compressedData[0] = ((unsigned char*)&shortWidth)[0];
 compressedData[1] = ((unsigned char*)&shortWidth)[1];
 compressedData[2] = ((unsigned char*)&shortHeight)[0];
 compressedData[3] = ((unsigned char*)&shortHeight)[1];
 compressedData[4] = charBpp;

 d = 5;

 for (s=0; s<sSize; s++){
   // read bytes in, counting the number of identical consecutive bytes

   // the difference
   diff = pixelData[s] - previous[s];

   // peek at the next bytes until they're different or we reach the end of
   // the data stream
   topt = s+sameByteCount;
   if (topt/*s+sameByteCount*/ < sSize){
     // move out of loop: topt = s+sameByteCount;
     //currentFrameTemp =  peek(s+sameByteCount, pixelData, sSize);
     //previousFrameTemp = peek(s+sameByteCount, previous, sSize);
     diff2 = pixelData[topt] - previous[topt];

     while (diff2 == diff && (sameByteCount < 254) && (topt/*s+sameByteCount*/+1 < sSize)){
       sameByteCount++;

       //topt = s+sameByteCount;
       topt++;
       //currentFrameTemp = peek(s+sameByteCount, pixelData, sSize);
       //previousFrameTemp = peek(s+sameByteCount, previous, sSize);
       diff2 = pixelData[topt] - previous[topt];
     }
   }

   // encode the 9 bit diff to 8 bits by throwing away lsb
   if (diff < 0){//pixelData[s] < previous[s]){
     // then we have a negative difference
     //diff = previous[s] - pixelData[s]; // get the positive difference
     diff = -diff; // (faster) get the positive difference

     //diff = diff & 255; // diff can at most be 255 - don't worry about top bits
     diff = diff >> 1;               // throw away lsb
     diff = diff | 128;              // put a sign bit in
   }
   else{
     // else we have a positive or zero difference
     //diff = diff & 255; // diff can at most be 255 - don't worry about top bits
     diff = diff >> 1;               // throw away lsb
                                     // don't bother with sign bit
   }

   // encode data to stream

   if (sameByteCount == 1){
     // check to see if the byte is the escape char
     if (diff != ESCAPE_CHAR){
       compressedData[d] = (unsigned char)diff;
       d++;
     }
     else {
       compressedData[d] = ESCAPE_CHAR;
       compressedData[d+1] = 1;
       compressedData[d+2] = (unsigned char)diff;
       d+=3;
     }
   }
   else if (sameByteCount > 2){ // stick this case in here for superscale opt
     compressedData[d] = ESCAPE_CHAR;
     compressedData[d+1] = sameByteCount;
     compressedData[d+2] = (unsigned char)diff;
     d+=3;
     s+=sameByteCount - 1;
     sameByteCount = 1;
   }
   else if (sameByteCount == 2 && diff != ESCAPE_CHAR){ // leave testing diff
     compressedData[d] = (unsigned char)diff;
     compressedData[d+1] = (unsigned char)diff;
     d+=2;
     s+=1; // (sameByteCount - 1) == (2-1) == 1
     sameByteCount = 1;
   }
   else { // catch the case where we have 2 bytes, both escape char
     compressedData[d] = ESCAPE_CHAR;
     compressedData[d+1] = sameByteCount;
     compressedData[d+2] = (unsigned char)diff;
     d+=3;
     s+=sameByteCount - 1;
     sameByteCount = 1;
   }

 }

 *compressedSize = d;

 return EVIZ_SUCCESS;
}

int eViz_RLE_readMemory(const unsigned char *compressedData,
                       const unsigned int compressedSize,
                       unsigned int *width, unsigned int *height,
                       unsigned int *bpp, unsigned char *pixelData)
{
 unsigned int s=0, d=0; // counters for the IO data
 unsigned int i;

 unsigned short shortWidth, shortHeight;
 unsigned char charBpp;

 unsigned char t;

 if (compressedData == NULL || width == NULL
     || height == NULL || bpp == NULL || pixelData == NULL)
 {
   fprintf(stderr, "eViz_RLE_readMemory: one or more arguments was NULL\n");
   return EVIZ_ERROR;
 }

 // firstly retrieve the image dimensions
 ((unsigned char*)&shortWidth)[0] = compressedData[0];
 ((unsigned char*)&shortWidth)[1] = compressedData[1];
 ((unsigned char*)&shortHeight)[0] = compressedData[2];
 ((unsigned char*)&shortHeight)[1] = compressedData[3];
 charBpp = compressedData[4];

 // then unencode the rest of the data
 for (s=5; s<compressedSize; s++){
   t = compressedData[s];

   if (t != ESCAPE_CHAR){
     pixelData[d] = t;
     d++;
   }
   else {
     for (i=0; i<compressedData[s+1]; i++){
       pixelData[d+i] = compressedData[s+2];
     }
     d+=compressedData[s+1];
     s+=2;
   }
 }

 *width = shortWidth;
 *height = shortHeight;
 *bpp = charBpp;

 return EVIZ_SUCCESS;
}

int eViz_RLED_readMemory(const unsigned char* compressedData,
                        const unsigned char* previous,
                        const unsigned int compressedSize,
                        unsigned int* width, unsigned int* height,
                        unsigned int* bpp, unsigned char* pixelData)
{
 unsigned int s=0, d=0; // counters for the IO data
 unsigned int i;

 unsigned short shortWidth, shortHeight;
 unsigned char charBpp;

 unsigned char t;
 unsigned short diff;
 int sign;

 if (compressedData == NULL || previous == NULL
     || width == NULL || height == NULL || bpp == NULL || pixelData == NULL)
 {
   fprintf(stderr, "eViz_RLED_readMemory: one or more arguments was NULL\n");
   return EVIZ_ERROR;
 }

 // firstly retrieve the image dimensions
 ((unsigned char*)&shortWidth)[0] = compressedData[0];
 ((unsigned char*)&shortWidth)[1] = compressedData[1];
 ((unsigned char*)&shortHeight)[0] = compressedData[2];
 ((unsigned char*)&shortHeight)[1] = compressedData[3];
 charBpp = compressedData[4];

 // then unencode the rest of the data
 for (s=5; s<compressedSize; s++){
   t = compressedData[s];

   if (t != ESCAPE_CHAR){
     sign = t & 128;
     t = t << 1;
     if (sign == 128){
       if (previous[d] >= t)
         pixelData[d] = previous[d] - t;
       else
         pixelData[d] = 0;
     }
     else{
       if (previous[d] + t <= 255)
         pixelData[d] = previous[d] + t;
       else
         pixelData[d] = 255;
     }
     d++;
   }
   else {
     diff = compressedData[s+2];
     sign = diff & 128;
     diff = diff << 1;
     diff = diff & 255;
     for (i=0; i<compressedData[s+1]; i++){
       if (sign == 128){
         if (previous[d+i] >= diff)
           pixelData[d+i] = previous[d+i] - diff;
         else
           pixelData[d+i] = 0;
       }
       else{
         if (previous[d+i] + diff <= 255)
           pixelData[d+i] = previous[d+i] + diff;
         else
           pixelData[d+i] = 255;
       }
     }
     d+=compressedData[s+1];
     s+=2;
   }
 }

 *width = shortWidth;
 *height = shortHeight;
 *bpp = charBpp;

 return EVIZ_SUCCESS;
}

/* ***************************************************************** */
/* Local function definitions                                        */
/* ***************************************************************** */

// roll this back in when everything works - its here for
// reasons of conceptual clarity for now
/* int peek(int i, unsigned char *array, int arraySize) */
/* { */
/*   if (i>=arraySize) */
/*     return -1; */
/*   else  */
/*     return array[i]; */
/* } */

// test.cpp : Defines the entry point for the console application.
//

