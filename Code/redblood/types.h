// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

