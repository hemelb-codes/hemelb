/*----------------------------------------------------------------------------

   (C)Copyright 2006 The University of Manchester, United Kingdom,
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

   Initial version by: M Riding, 05.04.2006

---------------------------------------------------------------------------*/

#ifndef __EVIZTYPES_H__
#define __EVIZTYPES_H__

#ifdef __cplusplus
extern "C" {
#endif

/* Function return values */
#define EVIZ_ERROR -1
#define EVIZ_SUCCESS 0

/* Big or little endian flags */
#ifdef __sgi
#define EVIZ_BIG_ENDIAN
#else
#define EVIZ_LITTLE_ENDIAN
#endif


/** Default port for images to be returned on. Can be overridden by setting
 *  the EVIZ_OUTPUT_STREAM_PORT environment variable.
 */
#define EVIZ_OUTPUT_PORT  50505

/** Default port for control data to be sent on. Can be overridden by setting
 *  the EVIZ_CONTROL_PORT environment varaiable.
 */
#define EVIZ_CONTROL_PORT 50506

/** Default port to communicate with the broker on. Can be overridden by
 *  setting the EVIZ_BROKER_PORT environment variable.
 */
#define EVIZ_BROKER_PORT  50507

/** Default port for gViz communication. Can be overridden by setting the
 *  EVIZ_GVIZ_PORT environment variable.
 */
#define EVIZ_GVIZ_PORT    50508


/** Output stream type: video (sequence of images)
 */
#define VIDEO_OUTPUT_STREAM 0

/** Output stream type: vrml (not implemented)
 */
#define VRML_OUTPUT_STREAM 1


/** Output stream protocol version: Major
 */
#define OUTPUT_STREAM_PROTOCOL_MAJOR_VERSION 1

/** Output stream protocol version: Minor
 */
#define OUTPUT_STREAM_PROTOCOL_MINOR_VERSION 8


/** Maximum size of data to read from a socket in one go
 */
#define MAX_DATA_SIZE 1024

/** Maximum number of connections on a socket
 */
#define EVIZ_BACKLOG 10

/* #define ARCH_LINUX_x86 */
/* #define ARCH_LINUX_x64 */
/* #define ARCH_LINUX_PPC */

/** Operating system: Linux
 */
#define ARCH_LINUX    0

/** Operating system: Win 32
 */
#define ARCH_WIN32    1

/** Operating sytem: MacOS
 */
#define ARCH_MAC      2

/** Operating system: HPUX
 */
#define ARCH_HPUX     3

/** Operating system: Irix
 */
#define ARCH_IRIX     4

/** Operating system: Solaris
 */
#define ARCH_SOLARIS  5

/** Operating system: AIX
 */
#define ARCH_AIX      6

/** Operating system: FreeBSD
 */
#define ARCH_FREEBSD  7

/** Operating system: Win 64
 */
#define ARCH_WIN64    8


/** Server software type: test_app
 */
#define TEST_APP      0

/** Server software type: AVS/Express
 */
#define AVS           10

/** Server software type: IRIS Explorer
 */
#define IRIS_EXPLORER 11

/** Server software type: OpenDX
 */
#define OPEN_DX       12

/** Server software type: RTRT
 */
#define RTRT          13

/** Server software type: VLib
 */
#define VLIB          14

/** Server software type: VTK
 */
#define VTK           15

/** Server software type: VMD
 */
#define VMD           16


/** Do not perform any image compression.
 */
#define NONE   0

/** Color cell compression codec. A fixed 8:1 compression codec.
 */
#define CCC    1

/** JPEG compression. Quality set to 75%.
 */
#define JPEG   2

/** PNG compression.
 */
#define PNG    3

/** Run time length encoding of difference images.
 */
#define RLED   4

/** Render locally.
 */
#define DIRECT 10


#define EVIZ_STRING        0
#define EVIZ_INTEGER       1
#define EVIZ_REAL          2
#define EVIZ_STRING_ARRAY  4
#define EVIZ_INTEGER_ARRAY 5
#define EVIZ_REAL_ARRAY    6
// Bit mask for testing array type
#define EVIZ_ARRAY         4

#define EVIZ_VIEW 0
#define EVIZ_STEER 1

typedef union eVizParameterValueUnion {
 char*          stringValue;
 int            integerValue;
 float          realValue;
 char**         stringArrayValue;
 int*           integerArrayValue;
 float*         realArrayValue;
} eVizParameterValueType;

typedef struct eVizParameterStruct {
 char name[64];
 int dataType; // string, integer, real, or array of string, integer etc.
 int steerType;
 int numberOfElements; // 1 for scalar > 1 for array/vector
 int hasMinimumValue;
 int hasMaximumValue;
 int nonQueueing;
 eVizParameterValueType minimumValue;
 eVizParameterValueType maximumValue;
 eVizParameterValueType currentValue;
} eVizParameterType;

#ifdef __cplusplus
}
#endif

#endif

