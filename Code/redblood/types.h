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

#include <set>
#include <utility>

#include "redblood/CellBase.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    //! Typical cell container type
    typedef std::set<std::shared_ptr<CellBase>> CellContainer;
    //! \brief Container of template meshes
    //! \details An instance of this object is used to reference meshes across the simulation
    typedef std::map<std::string, std::shared_ptr<CellBase>> TemplateCellContainer;
    //! Function to insert cells somewhere
    typedef std::function<void(CellContainer::value_type)> CellInserter;
  }
} // namespace hemelb::redblood
#endif

