// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_COMMON_INDEX_H
#define HLBGMYTOOL_COMMON_INDEX_H

#include <exception>
#include "util/Vector3D.h"

class IndexError : public std::exception {
 public:
  virtual const char* what() const noexcept { return "IndexError"; }
};

typedef hemelb::util::Vector3D<int> Index;
typedef hemelb::util::Vector3D<double> Vector;

#endif  // HLBGMYTOOL_COMMON_INDEX_H
