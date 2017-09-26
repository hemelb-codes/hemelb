// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELBSETUPTOOL_TRITREE_H
#define HEMELBSETUPTOOL_TRITREE_H

#include "Oct.h"
#include <boost/container/flat_set.hpp>

typedef boost::container::flat_set<int> IdList;
typedef Octree<IdList> TriTree;

std::ostream& operator<<(std::ostream& os, const IdList& lst);

#endif
