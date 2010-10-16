#ifndef HEMELB_DEBUG_PLATFORM_H
#define HEMELB_DEBUG_PLATFORM_H

/* This file hides the various platform differences.  Need to
 * cross-check definitions of the various macros against the Rules.mk
 * file.
 */

#ifdef HEMELB_DBG_LIBRARY_OSX
#include "dbg/OSX/OsxDebugger.h"

#elif defined(HEMELB_DBG_LIBRARY_LINUX)
#include "dbg/linux/LinuxDebugger.h"

#else
#include "dbg/none/NullDebugger.h"

#endif

#endif // HEMELB_DEBUG_PLATFORM_H
