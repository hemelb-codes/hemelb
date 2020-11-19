// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELBSETUPTOOL_INDEX_H
#define HEMELBSETUPTOOL_INDEX_H

#include <exception>
#include "util/Vector3D.h"

class IndexError: public std::exception {
	virtual const char* what() const throw () {
		return "IndexError";
	}
};

typedef hemelb::util::Vector3D<int> Index;
typedef hemelb::util::Vector3D<double> Vector;

#endif // HEMELBSETUPTOOL_INDEX_H
