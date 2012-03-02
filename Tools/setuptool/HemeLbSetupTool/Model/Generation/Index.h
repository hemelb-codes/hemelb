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
