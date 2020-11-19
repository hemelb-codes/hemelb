// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_PARMETISHEADER_H
#define HEMELB_GEOMETRY_PARMETISHEADER_H

#include "parmetis.h"
#if (PARMETIS_MAJOR_VERSION < 4)
typedef idxtype idx_t;
typedef float real_t;
#endif

#endif /* HEMELB_GEOMETRY_PARMETISHEADER_H */
