// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_DEBUG_PLATFORM_H
#define HEMELB_DEBUG_PLATFORM_H

/* This file hides the various platform differences.  Need to
 * cross-check definitions of the various macros against the Rules.mk
 * file.
 */

#ifdef HEMELB_DEBUG_LIBRARY_OSX
#include "debug/OSX/OsxDebugger.h"

#elif defined(HEMELB_DEBUG_LIBRARY_LINUX)
#include "debug/linux/LinuxDebugger.h"

#else
#include "debug/none/NullDebugger.h"

#endif

#endif // HEMELB_DEBUG_PLATFORM_H
