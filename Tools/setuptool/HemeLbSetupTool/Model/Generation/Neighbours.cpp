#include "Neighbours.h"

#include "io/formats/geometry.h"
// shortcut to geometry class
using hemelb::io::formats::geometry;

void Neighbours::Init() {
	// copy in the vectors of the neighbourhood
	Neighbours::vectors = geometry::GetNeighbourhood();

	//unsigned int iLater = 0;
	for (unsigned int i = 0, iLater = 0; i < Neighbours::n; ++i) {
		Neighbours::norms[i] = Neighbours::vectors[i].GetMagnitude();

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

		// Complete the table of inverse indices.
		for (unsigned int j = i; j < Neighbours::n; ++j) {
			if (Neighbours::vectors[i] == Neighbours::vectors[j]) {
				Neighbours::inverses[i] = j;
				Neighbours::inverses[j] = i;
				break;
			}
		}
	}
}

std::vector<unsigned int> Neighbours::laterNeighbourIndices(Neighbours::nLater);
std::vector<Index> Neighbours::vectors(Neighbours::n);
std::vector<unsigned int> Neighbours::inverses(Neighbours::n);
std::vector<double> Neighbours::norms(Neighbours::n);

