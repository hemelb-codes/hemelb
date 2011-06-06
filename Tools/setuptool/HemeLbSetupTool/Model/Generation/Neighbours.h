#ifndef HEMELBSETUPTOOL_NEIGHBOURS_H
#define HEMELBSETUPTOOL_NEIGHBOURS_H

#include <vector>

#include "D3Q15.h"
#include "Index.h"

struct Neighbours {
	// This assumes that the hemelb lattice descriptor class has a zero vector.

	enum {
		n = hemelb::D3Q15::NUMVECTORS - 1
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
