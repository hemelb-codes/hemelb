// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_HASSERT_H
#define HEMELB_HASSERT_H

#include "Exception.h"
#include "build_info.h"

// A HemeLB assertion macro
#define HASSERT(expr) \
if constexpr (::hemelb::DEBUG) { \
    if (!(expr)) { \
        throw (::hemelb::Exception() << "Assertion failure '" #expr "' in '" __FILE__ ":" << __LINE__); \
    } \
}

#endif
