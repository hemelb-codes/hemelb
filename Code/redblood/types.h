//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_TYPES_H
#define HEMELB_REDBLOOD_TYPES_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "redblood/types_fwd.h"
#include "redblood/CellBase.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    namespace details
    {
      //! Stable comparison of cells across nodes
      bool CellUUIDComparison::operator()(std::shared_ptr<CellBase> const&a,
					  std::shared_ptr<CellBase> const &b) const
      {
	return a->GetTag() < b->GetTag();
      }
    }
  }
} // namespace hemelb::redblood
#endif

