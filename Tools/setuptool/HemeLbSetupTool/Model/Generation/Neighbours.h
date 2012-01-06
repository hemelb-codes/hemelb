#ifndef HEMELBSETUPTOOL_NEIGHBOURS_H
#define HEMELBSETUPTOOL_NEIGHBOURS_H

#include <vector>

#include "lb/lattices/D3Q27.h"
#include "Index.h"

struct Neighbours {
	// This assumes that the hemelb lattice descriptor class has a zero vector.

	enum {
		n = hemelb::lb::lattices::D3Q27::NUMVECTORS - 1
	};
	enum {
		nLater = n / 2
	};
	static std::vector<unsigned int> laterNeighbourIndices;
	static std::vector<Index> vectors;
	static std::vector<double> norms;
	static std::vector<unsigned int> inverses;

	static void Init();
private:
	// Private to ensure it's not instantiated.
	Neighbours();
};
#endif // HEMELBSETUPTOOL_NEIGHBOURS_H
