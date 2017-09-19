// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_FLOODFILL_H
#define HEMELBSETUPTOOL_FLOODFILL_H

#include "FluidSiteTree.h"
#include "MaskTree.h"

// This function flood fills the domain
// And (eventually) propagates the fluid site count of children up the
// tree
class FloodFill {
public:
  typedef std::array<FluidTree::Int, 4> Idx;
  
  FloodFill(const FluidTree& t);
  
  MaskTree operator()() const;
  
  Idx GetStart() const;

private:
  
  class StopIteration : public std::exception {
  };

  const FluidTree& tree;
};
#endif
