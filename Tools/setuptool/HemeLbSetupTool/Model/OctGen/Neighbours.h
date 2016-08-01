// -*- mode: c++; -*-
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
