// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_FLOODFILL_H
#define HEMELBSETUPTOOL_FLOODFILL_H

#include "FluidSiteTree.h"

// This class flood fills the domain and propagates the fluid site count
// of children up the tree

void FloodFill(FluidTree& tree);

#endif
