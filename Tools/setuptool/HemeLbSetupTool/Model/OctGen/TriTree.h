// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_TRITREE_H
#define HEMELBSETUPTOOL_TRITREE_H

#include "Oct.h"
#include <boost/container/flat_set.hpp>

typedef boost::container::flat_set<int> IdList;
typedef Octree<IdList> TriTree;

#endif
