// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_FLOODFILL_H
#define HEMELBSETUPTOOL_FLOODFILL_H

#include "FluidSiteTree.h"

// This function flood fills the domain

// And (eventually) propagates the fluid site count of children up the
// tree

typedef Octree<bool> MaskTree;

MaskTree FloodFill(const FluidTree& tree);

typedef std::array<FluidTree::Int, 4> Idx;

class StopIteration : public std::exception {
};

Idx GetStart(const FluidTree& tree);

#endif
