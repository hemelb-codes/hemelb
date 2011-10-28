#include "Neighbours.h"

void Neighbours::Init() {
	// This assumes that the hemelb lattice descriptor class has a zero vector at the first entry.
	unsigned int iLater = 0;
	for (unsigned int i = 0; i < Neighbours::n; ++i) {
		Neighbours::vectors[i] = Index(hemelb::D3Q15::CX[i + 1],
				hemelb::D3Q15::CY[i + 1], hemelb::D3Q15::CZ[i + 1]);
		Neighbours::norms[i] = Neighbours::vectors[i].Magnitude<double>();

		// Figure out if this takes us to a later site
		if (Neighbours::vectors[i][0] > 0) {
			// yes
			Neighbours::laterNeighbourIndices[iLater] = i;
			iLater++;
		} else if (Neighbours::vectors[i][0] == 0) {
			// Maybe
			if (Neighbours::vectors[i][1] > 0) {
				// yes
				Neighbours::laterNeighbourIndices[iLater] = i;
				iLater++;
			} else if (Neighbours::vectors[i][1] == 0) {
				// Maybe
				if (Neighbours::vectors[i][2] > 0) {
					// yes
					Neighbours::laterNeighbourIndices[iLater] = i;
					iLater++;
				} else {
					// No
				}
			} else {
				// No
			}
		} else {
			// No
		}
		Neighbours::inverses[i] = hemelb::D3Q15::INVERSEDIRECTIONS[i + 1] - 1;
	}
}

std::vector<unsigned int> Neighbours::laterNeighbourIndices(Neighbours::nLater);
std::vector<Index> Neighbours::vectors(Neighbours::n);
std::vector<unsigned int> Neighbours::inverses(Neighbours::n);
std::vector<double> Neighbours::norms(Neighbours::n);

