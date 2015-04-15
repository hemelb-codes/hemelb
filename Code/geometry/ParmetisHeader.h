//
// Copyright (C) University College London, 2007-2013, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_GEOMETRY_PARMETISHEADER_H
#define HEMELB_GEOMETRY_PARMETISHEADER_H

#include "parmetis.h"
#if (PARMETIS_MAJOR_VERSION < 4)
typedef idxtype idx_t;
typedef float real_t;
#endif

#endif /* HEMELB_GEOMETRY_PARMETISHEADER_H */
