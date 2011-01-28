/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * par_metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: parmetislib.h,v 1.2 2003/07/21 17:50:22 karypis Exp $
 */

/*
#define DEBUG			1
#define DMALLOC			1
*/

#include "parmetis/stdheaders.h"
#include "parmetis/parmetis.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include "parmetis/rename.h"
#include "parmetis/defs.h"
#include "parmetis/struct.h"
#include "parmetis/macros.h"
#include "parmetis/proto.h"

