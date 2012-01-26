#ifndef HEMELBSETUPTOOL_NEIGHBOURS_H
#define HEMELBSETUPTOOL_NEIGHBOURS_H

#include <vector>
#include "Index.h"

#include "io/formats/geometry.h"
// shortcut to geometry class
using hemelb::io::formats::geometry;

struct Neighbours {
	// This assumes that the hemelb lattice descriptor class has a zero vector.

	enum {
		n = geometry::NumberOfDisplacements
	};
	enum {
		nLater = n / 2
	};

	static std::vector<unsigned int> laterNeighbourIndices;
	static geometry::DisplacementVector vectors;
	static std::vector<double> norms;
	static std::vector<unsigned int> inverses;

	static void Init();
private:
	// Private to ensure it's not instantiated.
	Neighbours();
};
#endif // HEMELBSETUPTOOL_NEIGHBOURS_H
