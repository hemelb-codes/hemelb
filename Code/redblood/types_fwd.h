//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_TYPES_FWD_H
#define HEMELB_REDBLOOD_TYPES_FWD_H

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>

namespace hemelb
{
  namespace redblood
  {
    class CellBase;

    namespace details
    {
      //! Stable comparison of cells across nodes
      struct CellUUIDComparison
      {
	inline bool operator()(std::shared_ptr<CellBase> const&a,
			std::shared_ptr<CellBase> const &b) const;
      };
    }
    //! Typical cell container type
    using CellContainer = std::set<std::shared_ptr<CellBase>, details::CellUUIDComparison>;
    //! \brief Container of template meshes
    //! \details An instance of this object is used to reference meshes across the simulation
    using TemplateCellContainer = std::map<std::string, std::shared_ptr<CellBase>> ;
    //! Function to insert cells somewhere
    using CellInserter = std::function<void(CellContainer::value_type)>;
  }
} // namespace hemelb::redblood
#endif

