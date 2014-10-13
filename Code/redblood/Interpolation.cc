//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/interpolation.h"

namespace hemelb { namespace redblood {

  void RegionIterator :: operator++() {
#     ifndef NDEBUG
        if(not isValid())
          throw Exception() << "Cannot increment invalid iterator\n";
#     endif
      if(current_[2] < max_[2]) {
        ++current_[2];
        return;
      }
      current_[2] = min_[2];
      if(current_[1] < max_[1]) {
        ++current_[1];
        return;
      }
      current_[1] = min_[1];
      ++current_[0];
  }

}}
