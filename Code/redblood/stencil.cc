//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/stencil.h"

namespace hemelb
{
  namespace redblood
  {
    namespace stencil
    {
#define HEMELB_STENCIL_MACRO(NAME, RANGE) const size_t NAME::range = RANGE;
      HEMELB_STENCIL_MACRO(FourPoint, 4u);
      HEMELB_STENCIL_MACRO(CosineApprox, 4u);
      HEMELB_STENCIL_MACRO(ThreePoint, 3u);
      HEMELB_STENCIL_MACRO(TwoPoint, 2u);
#undef HEMELB_STENCIL_MACRO
    }
  }
}  // hemelb::redblood::stencil
