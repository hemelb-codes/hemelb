// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELBSETUPTOOL_NEIGHBOURS_H
#define HEMELBSETUPTOOL_NEIGHBOURS_H

#include <vector>
#include "Vector.h"

class Neighbours {
public:
  static const Neighbours& Get();
  static const std::vector<Index>& GetDisplacements();
  static const std::vector<int>& GetOpposites();
  
  Neighbours(const Neighbours&) = delete;
  Neighbours& operator=(const Neighbours&) = delete;
  
private:
  Neighbours();
  std::vector<Index> neighbours;
  std::vector<int> inverses;
};
#endif
